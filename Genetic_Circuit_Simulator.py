
""" Authors:    Hamid Doosthosseini < hdoosth@mit.edu > , Voigt Lab , MIT
                Amin Espah Borujeni < amin.espah@gmail.com > , Voigt Lab , MIT """

""" Last updated: 06/25/2020"""

# ---------------------------------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------------------------------- #

import numpy as np
import math
import pandas as pd
import csv
import warnings
import scipy
from scipy import optimize
import pickle
from datetime import datetime
import regex as re
import string
import time
import os

gamma = 0.0067                                      # measured mRNA degradation rate
delta = np.log(2)/(45*60)                           # measured dilution rate
chemicals = []
dnacopy = {'Dv_pBM3R1-YFP':4.0,'0x41v70':9.0}   # measured DNA copy of plasmids

class chemical:
    def __init__(self, name, value=0.0, unit='binary'):
        self.name = name
        self.value = value
        self.unit=unit
        self.reset_value = value
    
    def set_value(self, value, unit):
        self.value = value
        self.unit=unit

    def reset(self):
        self.value = self.reset_value
    
    def make_current_reset(self):
        self.reset_value = self.value

class genetic_part:
    def __init__(self, genome='', part_type='null', start=0.0, end=0.0, strand='', info_dict={}):
        self.type = part_type
        self.strand = strand
        self.start = start
        self.end = end
        self.genome = genome
        self.repressed=False
        self.induced=False
        self.flux = 0.0
        self.added_flux = 0.0
        self.rna = 0.0
        self.protein = 0.0
        self.reset_flux = 0.0
        self.reset_added_flux = 0.0
        self.reset_rna = 0.0
        self.reset_protein = 0.0

        try:
            self.name = info_dict['Name']
        except:
            self.name = ''
            warnings.warn("Issue in resolving name for genetic part", UserWarning, stacklevel=2)
        if part_type == 'promoter':
            try:
                self.strength = float(info_dict['Strength'])
            except:
                self.strength = 0.0
                warnings.warn("Issue in resolving promoter parameters " + self.name, UserWarning, stacklevel=2)
            if 'Repressor' in info_dict.keys():
                try:
                    self.repressor_name = info_dict['Repressor']
                    self.repressed=True
                    self.K = max(float(info_dict['K']),1.0e-5)
                    self.n = float(info_dict['n'])
                    self.ymin = float(info_dict['ymin'])                
                except:
                    self.K = math.inf
                    self.n = 1.0
                    self.ymin = 1.0e-7
                    self.repressed=False
                    warnings.warn("Issue in resolving repressor parameters for " + self.name, UserWarning, stacklevel=2)
            if 'Inducer' in info_dict.keys():
                try:
                    self.K = max(float(info_dict['K']),1.0e-5)
                    self.ymin = float(info_dict['ymin'])
                    self.inducer_name = info_dict['Inducer']
                    self.induced=True
                except:
                    self.K = 0.0
                    self.ymin = 1.0e-7
                    self.induced=False
                    warnings.warn("Issue in resolving inducer parameters for " + self.name, UserWarning, stacklevel=2)

        elif part_type == 'terminator':
            try:
                self.strength = float(info_dict['Strength'])
            except:
                self.strength = 0.0
                warnings.warn("Issue in resolving terminator strength for " + self.name, UserWarning, stacklevel=2)

        elif part_type == 'ribozyme':
            try:
                self.efficiency = float(info_dict['Cleavage_Efficiency'])
                self.cut_site = int(info_dict['cut_site'])
            except:
                self.efficiency = 0.0
                warnings.warn("Issue in resolving ribozyme parameters for " + self.name, UserWarning, stacklevel=2)
        
        elif part_type == 'RBS':
            try:
                self.efficiency = float(info_dict['Translation_Efficiency'])
            except:
                self.efficiency = 0.0
                warnings.warn("Issue in resolving RBS translation efficiency for " + self.name, UserWarning, stacklevel=2)

    def set_neighbours(self, us, ds):
        self.us = us
        self.ds = ds

    def calculate(self):
        if self.type == 'gene':
            self.new_flux = self.us.flux
        elif self.type == 'ribozyme':
            self.new_flux = self.us.flux
        elif self.type == 'RBS':
            self.new_flux = self.us.flux
        elif self.type == 'terminator':
            self.new_flux = self.us.flux/self.strength
        elif self.type == 'promoter':
            
            if self.induced:
                self.added_flux = self.ymin + (self.strength-self.ymin)*self.inducer.value
                self.new_flux = self.us.flux + self.added_flux
            if self.repressed:
                self.added_flux = self.ymin + (self.strength-self.ymin)/(1+(sum([x.protein for x in self.repressors])/self.K)**self.n)
                self.new_flux = self.us.flux + self.added_flux
        else:
            self.new_flux = 0.0
    
    def update(self, dt):
        self.flux = self.new_flux
        if self.us.type == 'ribozyme':
            self.rna += dnacopy[self.genome]*self.flux*dt - (self.us.efficiency+2*(1-self.us.efficiency))*self.rna*gamma*dt
        else:
            self.rna += dnacopy[self.genome]*self.flux*dt - self.rna*gamma*dt
        if self.us.type == 'RBS':
            self.protein += self.us.rna*self.us.efficiency*dt - delta*self.protein*dt

    def reset(self):
        self.flux = self.reset_flux
        self.added_flux = self.reset_added_flux
        self.rna = self.reset_rna
        self.protein = self.reset_protein

    def make_current_reset(self):
        self.reset_flux = self.flux
        self.reset_rna = self.rna
        self.reset_protein = self.protein
        self.reset_added_flux = self.added_flux

def pickleloader(pklFile):
    try:
        while True:
            yield pickle.load(pklFile)
    except EOFError:
        pass

def parse_gff3(file_name):
    all_parts = []
    genomes = []
    with open(file_name) as file:
        readcsv = csv.reader(file)
        for row in readcsv:
            info = row[8].split(';')
            info_dict = {detail.split('=')[0]:detail.split('=')[1] for detail in info}
            part = genetic_part(row[0],row[2],float(row[3]),float(row[4]),row[6],info_dict)
            if row[6] == '+' or row[6] == '-':
                all_parts.append(part)
            else:
                warnings.warn("Warning: No parseable strand in input " + str(row), UserWarning, stacklevel=2)
            if row[0] not in genomes:
                genomes.append(row[0])
    
    for genome in genomes:
        all_parts.append(genetic_part(genome=genome,strand='+'))
        all_parts.append(genetic_part(genome=genome,strand='-'))
    all_parts = sorted(all_parts,key=lambda e: (e.genome,e.strand,e.start))
    all_num = len(all_parts)
    for i in range(all_num):
        all_parts[i].set_neighbours(all_parts[(i-1)%all_num],all_parts[(i+1)%all_num])
    for part in all_parts:
        if part.repressed:
            flag = True
            part.repressors = []
            for check_part in all_parts:
                if check_part.name == part.repressor_name:
                    part.repressors.append(check_part)
                    flag = False
            if flag:
                warnings.warn("Repressor " + str(part.repressor_name) + " (for " + str(part.name) + ") not found in circuit", UserWarning, stacklevel=2)
                part.repressed = False
    for part in all_parts:
        if part.induced:
            flag = True
            for inducer in chemicals:
                if inducer.name == part.inducer_name:
                    flag = False
                    part.inducer = inducer
            if flag:
                inducer = chemical(part.inducer_name)
                chemicals.append(inducer)
                part.inducer = inducer
                flag = False
            
    print('Loading circuit from GFF3 format complete: \n')
    print('List of parts in order including spacer entities: ' + str([part.name for part in all_parts]) + '\n')
    print('List of chemicals: ' + str([inducer.name for inducer in chemicals]) + '\n')
    return all_parts

def run_sim(all_parts,induction_state={},publish=False,runtime=12):

    for key,value in induction_state.items():
        count_inducers = 0
        for inducer in chemicals:
            if inducer.name == key:
                inducer.value = value
                count_inducers += 1
        if count_inducers < 1: 
            warnings.warn("No inducer " + key + " found in list of chemicals", UserWarning, stacklevel=2)
        if count_inducers > 1: 
            warnings.warn("Multiple entries for " + key + " found in list of chemicals", UserWarning, stacklevel=2)

    t = 0.0
    dt = 60
    data = {}
    for part in all_parts:
        if part.type == 'gene':
            data[part.name+'_mrna'] = [part.reset_rna,part.reset_rna]
            data[part.name+'_protein'] = [part.reset_protein,part.reset_protein]
        if part.type == 'promoter':
            data[part.name+'_flux'] = [part.reset_flux,part.reset_flux]
            data[part.name+'_addedflux'] = [part.reset_added_flux,part.reset_added_flux]
        if part.type == 'terminator':
            data[part.name+'_flux'] = [part.reset_flux,part.reset_flux]
        data['time'] = [-(60*60),-1e-5]
    for inducer in chemicals:
        data[inducer.name] = [inducer.reset_value,inducer.reset_value]
    
    while t<60*60*runtime:
        t += dt
        for part in all_parts:
            part.calculate()
        for part in all_parts:
            part.update(dt)
        data['time'].append(t)
        for part in all_parts:
            if part.type == 'gene':
                data[part.name+'_mrna'].append(part.rna)
                data[part.name+'_protein'].append(part.protein)
            if part.type == 'promoter':
                data[part.name+'_addedflux'].append(part.added_flux)
                data[part.name+'_flux'].append(part.flux)
            if part.type == 'terminator':
                data[part.name+'_flux'].append(part.flux)
        for inducer in chemicals:
            data[inducer.name].append(inducer.value)
    output_file = 'output/dynamic_run'+datetime.now().strftime("%Y_%m_%d_%H_%M_%S")+'.csv'
    if publish:
        pd.DataFrame.from_dict(data).to_csv(output_file)
        time.sleep(1)
    data_ss = {}
    for key in data:
        data_ss[key] = [data[key][-1]]
    return(data_ss,all_parts,output_file)

def sensitivity_analysis(fromfile=False):

    def circuit_score(all_parts,on_states,off_states):
        on_values = []
        off_values = []
        all_dat = []
        for induction_state in on_states:
            induction_state0={'IPTG':0.0, 'aTc':0.0, 'Ara':0.0}
            data0, all_parts, run_file = run_sim(all_parts,induction_state0)
            data, all_parts, run_file = run_sim(all_parts,induction_state)
            on_values.append(data['YFP_protein'][-1])
            all_dat.append({'induction_state':induction_state,'data':data})
            for inducer in chemicals:
                inducer.reset()
            for part in all_parts:
                part.reset()
        for induction_state in off_states:
            induction_state0={'IPTG':0.0, 'aTc':0.0, 'Ara':0.0}
            data0, all_parts, run_file = run_sim(all_parts,induction_state0)
            data, all_parts, run_file = run_sim(all_parts,induction_state)
            off_values.append(data['YFP_protein'][-1])
            all_dat.append({'induction_state':induction_state,'data':data})
            for inducer in chemicals:
                inducer.reset()
            for part in all_parts:
                part.reset()
        return max(min(on_values),1.0)/max(max(off_values),1.0), on_values, off_values, all_dat

    def sensitivity_analysis_ind(part,attr='param',fromfile=True,limits=[],x_s=[]):
        
        if fromfile:
            cur_value = getattr(part,attr)
            print(part.name, attr, cur_value, 'from file: ', fromfile)
            points = []
            for filename in os.listdir("./output/"):
                if filename.startswith("sa_"+part.name + "_" + attr) and filename.endswith('.dat'):
                    # print(filename)
                    with open("./output/"+filename,'rb') as handle:
                        for line in pickleloader(handle):
                            [name,x,y,on_values,off_values]  = line
                            if x>limits[0] and x<limits[1]:
                                # y_new = min(on_values)/max(max(off_values),dnacopy['Dv_pBM3R1-YFP']*1e-5/(gamma)*1.29/delta)
                                y_new = min(on_values)/max(off_values)
                                # points.append((x,y_new,on_values[0]/max(max(off_values),1.0/gamma),on_values[1]/max(max(off_values),1.0/gamma)))
                                points.append((x,y_new))
            points = sorted(points,key=lambda x: x[0])
            x_s = np.array([point[0] for point in points])
            y_s = [np.array([point[i+1] for point in points]) for i in range(len(points[0])-1)]
        
        else:
            cur_value = getattr(part,attr)
            print(part.name, attr, cur_value, 'from file: ', fromfile)
            if len(limits) == 0:
                limits = [cur_value*0.01,cur_value*100]
            if len(x_s) == 0:
                x_s = np.logspace(np.log10(limits[0]),np.log10(limits[1]),100)
                # x_s = np.linspace(limits[0],limits[1],1000)
            y_s = [[]]
            file_name = "output/sa_" + part.name + "_" + attr  + ".dat"
            with open(file_name,'ab') as handle:
                for x in x_s:
                    setattr(part,attr,x)
                    y, on_values, off_values, all_dat = circuit_score(all_parts,on_states,off_states)
                    y_s[0].append(y)
                    pickle.dump([part.name,x,y,on_values,off_values],handle)
                setattr(part,attr,cur_value)
        return(x_s,y_s)
    
    on_states = [{'IPTG':1,'aTc':0,'Ara':0},{'IPTG':1,'aTc':1,'Ara':1}]
    off_states = [{'IPTG':0,'aTc':0,'Ara':0},{'IPTG':0,'aTc':1,'Ara':0},
    {'IPTG':0,'aTc':0,'Ara':1},{'IPTG':1,'aTc':1,'Ara':0},
    {'IPTG':1,'aTc':0,'Ara':1},{'IPTG':0,'aTc':1,'Ara':1}]
    
    all_parts = parse_gff3('Circuit_Parametrized.csv')
    dict_sa = {}

    print("Circuit of current parameters : ", circuit_score(all_parts,on_states,off_states))
    for part in all_parts:
        dict_sa[part] = {}
        if part.type == 'promoter':
            x_s, y_s = sensitivity_analysis_ind(part,attr='strength',limits=[1e-5,1e2],fromfile=fromfile)
            dict_sa[part]['strength'] = (x_s, y_s)
            x_s, y_s = sensitivity_analysis_ind(part,attr='K',limits=[1e-3,1e7],fromfile=fromfile)#,limits=[1e-5,1e3]
            dict_sa[part]['K'] = (x_s, y_s)
            x_s, y_s = sensitivity_analysis_ind(part,attr='ymin',limits=[1.0e-7,1.0e-3],fromfile=fromfile)#,limits=[1e-5,1e3]
            dict_sa[part]['ymin'] = (x_s, y_s)
        if part.type == 'terminator':
            x_s, y_s = sensitivity_analysis_ind(part,attr='strength',limits=[1,1e4],fromfile=fromfile)
            dict_sa[part]['strength'] = (x_s, y_s)
        if part.type == 'ribozyme':
            x_s, y_s = sensitivity_analysis_ind(part,attr='efficiency',limits=[0,1],x_s=np.linspace(0,1,1000),fromfile=fromfile)
            dict_sa[part]['efficiency'] = (x_s, y_s)
        if part.type == 'RBS':
            x_s, y_s = sensitivity_analysis_ind(part,attr='efficiency',limits=[1e-3,1e3],fromfile=fromfile)
            dict_sa[part]['efficiency'] = (x_s, y_s)
    
    return dict_sa

def kn_fitting():
    
    x0 = [0.001,2.0]           # intial guesses
    k_calc_filename = 'gate_data.csv'
    K_data = pd.read_csv(k_calc_filename)
    print("\n-- Fitting K and n values from file: {} --".format(k_calc_filename))
    print("-- Initial guess for all gates: k [RNAP/s] = {:.2e}, n = {:.2f} --".format(x0[0],x0[1]))
   
    def K_err(variables,in_data,out_data,te):
        k,n = variables
        ymax = max(out_data)
        ymin = min(out_data)
        y_i = ymax - out_data
        x_i = np.log((in_data/gamma)*(te/delta))
        a_i = out_data-ymin
        beta = n
        alpha = -n*np.log((k/gamma)*(te/delta))
        return sum((y_i-a_i*np.exp(beta*x_i+alpha))*x_i)**2
                
    with open('kn_regression.txt', 'w') as handle:
        for col in K_data:
            if col.endswith('input'):
                in_data0 = K_data[col].values
                out_data0 = K_data[col.split('_')[0]+"_output"].values
                in_data = in_data0[(np.isfinite(in_data0)) & (np.isfinite(out_data0))]
                out_data = out_data0[(np.isfinite(in_data0)) & (np.isfinite(out_data0))]
                te = K_data[col.split('_')[0]+"_TE"].values[0]
                results = scipy.optimize.minimize(lambda x: K_err(x,in_data,out_data,te),x0,method='Nelder-Mead',tol=1e-12,options={'maxiter':1e5})
                print(col,results['x'][0]/0.0067*te/delta, results['x'][1], results['success'], results['fun'])
                print(col, results, file=handle)

# ---------------------------------- #

kn_fitting()            # Fitting of K, n to flux data. Input from 'gate_data.csv' file and output to 'kn_regression.txt'

# ---------------------------------- #	

# run_sim runs the simulation from the current situation of all genetic parts
# outputs results to csv file: './output/dynamic_run'+datetime.now().strftime("%Y_%m_%d_%H_%M_%S")+'.csv'
# Example: setting initial conditions to only IPTG induction and simulate induction with all 3

all_parts = parse_gff3('Circuit_Parametrized.csv')                # reads input GFF3 file
induction_state0={'IPTG':0.0, 'aTc':0.0, 'Ara':0.0}     # initial induction state
data0, all_parts, run_file = run_sim(all_parts,induction_state0,publish=True)   # simulates state with induction_state0
# data0: full data, all_parts: updated state of all genetic parts in the cell, run_file: name of file containing output data
 
for part in all_parts:                                  # set all parts and inductions to current steady state
    part.make_current_reset()
    part.reset()
for inducer in chemicals:
    inducer.make_current_reset()
    inducer.reset()

# Simulate with induction_state1 from prior steady state at induction_state0
induction_state1={'IPTG':1.0, 'aTc':1.0, 'Ara':1.0}
data_ss, all_parts, run_file = run_sim(all_parts,induction_state1,publish=True,runtime=12)

# ---------------------------------- #	

# conducts complete sensitivity and outputs to file in ./output/ in binary format. Returns dictionary of two lists for each part:
# x_s is the list of varied part parameter and y_s is the list of corresponding circuit scores
# calling the function with fromfile=True opens said files and returns the same dictionary

dict_sa = sensitivity_analysis(fromfile=False)  

# ---------------------------------- #	