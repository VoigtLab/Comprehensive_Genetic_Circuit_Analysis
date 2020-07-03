

""" Authors: 	Amin Espah Borujeni < amin.espah@gmail.com > , Voigt Lab , MIT
		Jing Zhang < jgzhang@bu.edu > , Voigt Lab , MIT """

""" Last updated: 06/25/2020"""

# ---------------------------------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------------------------------- #

import numpy as np
import matplotlib.pyplot as plt
import pickle
import copy
import xlrd
import xlwt
import math
from scipy import stats, optimize
import subprocess

# ---------------------------------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------------------------------- #

class Genetic_Circuit_Analysis_using_Sequening_Methods():
	
	
	def __init__(self, par, inducer, degradation, RNAP_p15A_total, RNAP_p15A, coef1, coef2, omega, total_active_ribosomes, codon_table, amino_acid_MW, frag_dist_bin, \
				reference_name, reference_length, Linker_Sequence, RAW_FASTQ_filename, Trimmed_RAW_FASTQ_filename, Unaligned_FASTQ_filename, \
				Coding_Genes_Fasta_filename, Coding_Genes_INDEXED_Fasta_Filename, Coding_Genes_GFF_filename, Coding_Genes_MAP_output_filename, \
				circuit_output_RD_exponential_fit_file, reporter_output_RD_exponential_fit_file, genome_output_RD_exponential_fit_file, tRNA_rRNA_operon_GFF_filename, \
				output_RNASeq_reads_filename, output_RiboSeq_reads_filename, output_RNA_seq_Profiles_filename, output_Ribo_seq_Profiles_filename, \
				output_FPKM_table_filename, output_RD_table_filename, output_Protein_MW_table_filename, output_Proteome_Fraction_table_filename):
				 
		self.par = par
		self.inducer = inducer
		self.degradation = degradation
		self.RNAP_p15A_total = RNAP_p15A_total
		self.RNAP_p15A = RNAP_p15A
		self.coef1 = coef1
		self.coef2 = coef2
		self.omega = omega
		self.total_active_ribosomes = total_active_ribosomes
		self.codon_table = codon_table
		self.amino_acid_MW = amino_acid_MW
		self.frag_dist_bin = frag_dist_bin
		
		self.reference_name = reference_name
		self.reference_length = reference_length
		self.Linker_Sequence = Linker_Sequence
		self.RAW_FASTQ_filename = RAW_FASTQ_filename
		self.Trimmed_RAW_FASTQ_filename = Trimmed_RAW_FASTQ_filename
		self.Unaligned_FASTQ_filename = Unaligned_FASTQ_filename
		self.Coding_Genes_Fasta_filename = Coding_Genes_Fasta_filename
		self.Coding_Genes_INDEXED_Fasta_Filename = Coding_Genes_INDEXED_Fasta_Filename
		self.Coding_Genes_GFF_filename = Coding_Genes_GFF_filename
		self.Coding_Genes_MAP_output_filename = Coding_Genes_MAP_output_filename
		
		self.circuit_output_RD_exponential_fit_file = circuit_output_RD_exponential_fit_file
		self.reporter_output_RD_exponential_fit_file = reporter_output_RD_exponential_fit_file
		self.genome_output_RD_exponential_fit_file = genome_output_RD_exponential_fit_file
		self.tRNA_rRNA_operon_GFF_filename = tRNA_rRNA_operon_GFF_filename
		
		self.output_RNASeq_reads_filename = output_RNASeq_reads_filename
		self.output_RiboSeq_reads_filename = output_RiboSeq_reads_filename
		self.output_RNA_seq_Profiles_filename = output_RNA_seq_Profiles_filename
		self.output_Ribo_seq_Profiles_filename = output_Ribo_seq_Profiles_filename
		
		self.output_FPKM_table_filename = output_FPKM_table_filename
		self.output_RD_table_filename = output_RD_table_filename
		self.output_Protein_MW_table_filename = output_Protein_MW_table_filename
		self.output_Proteome_Fraction_table_filename = output_Proteome_Fraction_table_filename
		
	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def reverse_complement_naive(self, seq):
		
		this_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
		
		return "".join(this_complement.get(base.upper(), base) for base in reversed(seq))

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def read_gff(self, inputfiles, item):
		
		gene_dict = {aa:{} for aa in item}
		
		for inputfile in inputfiles:
			
			lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
			
			for line in lines:
				tokens = line.split('\t')
				genome, feature, start, end, strand = tokens[0], tokens[2], int(tokens[3]), int(tokens[4]), tokens[6]
				name = (tokens[-1].split('Name='))[-1]
				gene_dict[genome][name] = {'start': start, 'end': end, 'strand': strand, 'feature': feature}
				
		return gene_dict

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def amino_acid_seq (self, seq):
	
		codon = [seq[i:i+3].upper() for i in range(0, len(seq), 3)]
		
		aminoacid_seq = []
		if codon[0] in ['ATG','GTG','CTG','TTG']:
			for item in codon:
				for key in self.codon_table:
					if item in self.codon_table[key]:
						aminoacid_seq.append(key)
						break

		return aminoacid_seq

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def model_func(self, t, A, K, C):
		
		return A * np.exp(-K * t) + C

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def map_raw_FASTQ_files(self):
	
		for sample_name in self.RAW_FASTQ_filename:
			
			#Linker removal
			#cutadapt (download site: https://cutadapt.readthedocs.org/en/stable/)
			cmd_index = 'cutadapt -a ' + self.Linker_Sequence + ' —-discard-untrimmed -m 10 -o ' + self.Trimmed_RAW_FASTQ_filename[sample_name] + ' ' + self.RAW_FASTQ_filename[sample_name]
			print("Removing linker: "+cmd_index)
			status0 = subprocess.call(cmd_index, shell=True)
			
			# ----------------------------------- #
			
			# Make the indexes
			cmd_index = 'bowtie-build' + ' -p ' + self.Coding_Genes_Fasta_filename[sample_name] + ' ' + self.Coding_Genes_INDEXED_Fasta_Filename[sample_name]
			print("Making index coding genes: "+cmd_index)
			status2 = subprocess.call(cmd_index, shell=True)
			
			# ----------------------------------- #
			
			# Perform the mapping
			cmd_mapping = 'bowtie' + ' —t -v1 -m2 -k1 —-un ' + self.Unaligned_FASTQ_filename[sample_name] + ' ' + self.Coding_Genes_INDEXED_Fasta_Filename[sample_name] + ' ' + self.Trimmed_RAW_FASTQ_filename[sample_name] + ' ' + self.Coding_Genes_MAP_output_filename[sample_name]
			print("Mapping Reads Bowtie: "+cmd_mapping)
			status4 = subprocess.call(cmd_mapping, shell=True)

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def generate_mapped_reads(self):
		
		RNA_seq_Read = {}
		Ribo_seq_Read = {}
		
		for sample_name in self.Coding_Genes_MAP_output_filename:
			
			if sample_name.split('_')[0] == 'RNASeq':
				Output_directory = self.output_RNASeq_reads_filename
				
			elif sample_name.split('_')[0] == 'RiboSeq':
				Output_directory = self.output_RiboSeq_reads_filename
			
			(circuit_plasmid_name, reporter_plasmid_name, genome_name) = self.reference_name[sample_name]
			(circuit_plasmid_length, reporter_plasmid_length, genome_length) = self.reference_length[sample_name]
			
			circuit_fwd_avg_read = {}
			circuit_rev_avg_read = {}
			for j in range(circuit_plasmid_length):
				circuit_fwd_avg_read[j] = [1]
				circuit_rev_avg_read[j] = [1]

			reporter_fwd_avg_read = {}
			reporter_rev_avg_read = {}
			for k in range(reporter_plasmid_length):
				reporter_fwd_avg_read[k] = [1]
				reporter_rev_avg_read[k] = [1]

			genome_fwd_avg_read = {}
			genome_rev_avg_read = {}
			for k in range(genome_length):
				genome_fwd_avg_read[k] = [1]
				genome_rev_avg_read[k] = [1]
			
			
			f_in = open(self.Coding_Genes_MAP_output_filename[sample_name], 'rU')
			lines = f_in.readlines()
			f_in.close()
			for line in lines:
				tokens = line.replace('\n', '').split('\t')
				tokens = filter(None, tokens)
				col2   = str(tokens[1])   # strand 
				col3   = str(tokens[2])   # Chromosome 
				col4   = int(tokens[3])   # left-most position
				
				try: 
					col5 = str(tokens[4])   # footprint seq 
				except IndexError: 
					print(line)
					
				length = len(col5)
				
				if col3 == circuit_plasmid_name:
				
					if col2 == '+':       # for plus strand
						if len(tokens) > 7: 
							if tokens[7].startswith('0'):   # if there is mismatch in the 1st position 
								length0 = length - 1        # substract wrong base at 1st position 
								end5 = col4 + 1             # substract 1st base
								end3 = end5 + length0 - 1
								circuit_fwd_avg_read[end5].append(length0)
							else: 
								end5 = col4
								end3 = end5 + length - 1
								circuit_fwd_avg_read[end5].append(length)
						else: 
							end5 = col4
							end3 = end5 + length - 1
							circuit_fwd_avg_read[end5].append(length)
							
					elif col2 == '-':
						if len(tokens) > 7:
							if tokens[7].startswith('0'):
								length0 = length - 1
								end3 = col4      # for minus strand, bowtie gives leftmost position (3' end) with zero-based numbering
								end5 = end3 + length0 - 1
								circuit_rev_avg_read[end5].append(length0)
							else: 
								end3 = col4
								end5 = end3 + length - 1
								circuit_rev_avg_read[end5].append(length)
						else: 
							end3 = col4
							end5 = end3 + length - 1
							circuit_rev_avg_read[end5].append(length)
				
				elif col3 == reporter_plasmid_name:
				
					if col2 == '+':       # for plus strand
						if len(tokens) > 7: 
							if tokens[7].startswith('0'):   # if there is mismatch in the 1st position 
								length0 = length - 1        # substract wrong base at 1st position 
								end5 = col4 + 1             # substract 1st base
								end3 = end5 + length0 - 1
								reporter_fwd_avg_read[end5].append(length0)
							else: 
								end5 = col4
								end3 = end5 + length - 1
								reporter_fwd_avg_read[end5].append(length)
						else: 
							end5 = col4
							end3 = end5 + length - 1
							reporter_fwd_avg_read[end5].append(length)
							
					elif col2 == '-':
						if len(tokens) > 7:
							if tokens[7].startswith('0'):
								length0 = length - 1
								end3 = col4      # for minus strand, bowtie gives leftmost position (3' end) with zero-based numbering
								end5 = end3 + length0 - 1
								reporter_rev_avg_read[end5].append(length0)
							else: 
								end3 = col4
								end5 = end3 + length - 1
								reporter_rev_avg_read[end5].append(length)
						else: 
							end3 = col4
							end5 = end3 + length - 1
							reporter_rev_avg_read[end5].append(length)
			
				elif col3 == genome_name:
				
					if col2 == '+':       # for plus strand
						if len(tokens) > 7: 
							if tokens[7].startswith('0'):   # if there is mismatch in the 1st position 
								length0 = length - 1        # substract wrong base at 1st position 
								end5 = col4 + 1             # substract 1st base
								end3 = end5 + length0 - 1
								genome_fwd_avg_read[end5].append(length0)
							else: 
								end5 = col4
								end3 = end5 + length - 1
								genome_fwd_avg_read[end5].append(length)
						else: 
							end5 = col4
							end3 = end5 + length - 1
							genome_fwd_avg_read[end5].append(length)
							
					elif col2 == '-':
						if len(tokens) > 7:
							if tokens[7].startswith('0'):
								length0 = length - 1
								end3 = col4      # for minus strand, bowtie gives leftmost position (3' end) with zero-based numbering
								end5 = end3 + length0 - 1
								genome_rev_avg_read[end5].append(length0)
							else: 
								end3 = col4
								end5 = end3 + length - 1
								genome_rev_avg_read[end5].append(length)
						else: 
							end3 = col4
							end5 = end3 + length - 1
							genome_rev_avg_read[end5].append(length)
							
			
			with open(Output_directory + '/circuit_fwd_avg_read__' + sample_name + '.txt', 'wb') as handle:
				pickle.dump(circuit_fwd_avg_read, handle)

			with open(Output_directory + '/circuit_rev_avg_read__' + sample_name + '.txt', 'wb') as handle1:
				pickle.dump(circuit_rev_avg_read, handle1) 

			with open(Output_directory + '/reporter_fwd_avg_read__' + sample_name + '.txt', 'wb') as handle2:
				pickle.dump(reporter_fwd_avg_read, handle2)

			with open(Output_directory + '/reporter_rev_avg_read__' + sample_name + '.txt', 'wb') as handle3:
				pickle.dump(reporter_rev_avg_read, handle3)

			with open(Output_directory + '/genome_fwd_avg_read__' + sample_name + '.txt', 'wb') as handle4:
				pickle.dump(genome_fwd_avg_read, handle4)

			with open(Output_directory + '/genome_rev_avg_read__' + sample_name + '.txt', 'wb') as handle5:
				pickle.dump(genome_rev_avg_read, handle5) 


			if sample_name.split('_')[0] == 'RNASeq':
				
				RNA_seq_Read[sample_name] = {circuit_plasmid_name:{},reporter_plasmid_name:{},genome_name:{}}
				
				RNA_seq_Read[sample_name][circuit_plasmid_name]['fwd'] = circuit_fwd_avg_read
				RNA_seq_Read[sample_name][circuit_plasmid_name]['rev'] = circuit_rev_avg_read
				
				RNA_seq_Read[sample_name][reporter_plasmid_name]['fwd'] = reporter_fwd_avg_read
				RNA_seq_Read[sample_name][reporter_plasmid_name]['rev'] = reporter_rev_avg_read
				
				RNA_seq_Read[sample_name][genome_name]['fwd'] = genome_fwd_avg_read
				RNA_seq_Read[sample_name][genome_name]['rev'] = genome_rev_avg_read
				
				
			elif sample_name.split('_')[0] == 'RiboSeq':
				
				Ribo_seq_Read[sample_name] = {circuit_plasmid_name:{},reporter_plasmid_name:{},genome_name:{}}
				
				Ribo_seq_Read[sample_name][circuit_plasmid_name]['fwd'] = circuit_fwd_avg_read
				Ribo_seq_Read[sample_name][circuit_plasmid_name]['rev'] = circuit_rev_avg_read
				
				Ribo_seq_Read[sample_name][reporter_plasmid_name]['fwd'] = reporter_fwd_avg_read
				Ribo_seq_Read[sample_name][reporter_plasmid_name]['rev'] = reporter_rev_avg_read
				
				Ribo_seq_Read[sample_name][genome_name]['fwd'] = genome_fwd_avg_read
				Ribo_seq_Read[sample_name][genome_name]['rev'] = genome_rev_avg_read
				
				
		return RNA_seq_Read, Ribo_seq_Read
	
	# ---------------------------------------------------------------------------------------------------------------------------------------- #
	
	def generate_RNA_Seq_profile(self):
		
		RNA_seq_total_read = {}
		RNA_seq_total_read_tRNA_rRNA_operon_deleted = {}
		RNA_seq_total_tRNA_fraction = {}
		
		RNA_seq_Profile = {}
		
		for sample_name in self.Coding_Genes_MAP_output_filename:
			
			if sample_name.split('_')[0] == 'RNASeq':
				
				Input_directory = self.output_RNASeq_reads_filename
			
				(circuit_plasmid_name, reporter_plasmid_name, genome_name) = self.reference_name[sample_name]
				(circuit_plasmid_length, reporter_plasmid_length, genome_length) = self.reference_length[sample_name]
				
				Master_wb_read = xlrd.open_workbook(self.tRNA_rRNA_operon_GFF_filename + '/NC_010473_1_operon.xls',"rb")
				Master_sh = Master_wb_read.sheet_by_name('NC_010473.1.operon')
				
				tRAN_rRNA_operon_range = []
				for i in range(43):
					left_bd = int(Master_sh.cell(2*i,11).value)
					right_bd = int(Master_sh.cell(2*i,12).value)
					tRAN_rRNA_operon_range.append((left_bd,right_bd))
				
				# ----------------------------------- #
				
				with open(Input_directory + '/circuit_fwd_avg_read__' + sample_name + '.txt', 'rb') as handle:
					circuit_fwd_avg_read_RNA = pickle.load(handle)

				with open(Input_directory + '/circuit_rev_avg_read__' + sample_name + '.txt', 'rb') as handle1:
					circuit_rev_avg_read_RNA = pickle.load(handle1)

				with open(Input_directory + '/reporter_fwd_avg_read__' + sample_name + '.txt', 'rb') as handle2:
					reporter_fwd_avg_read_RNA = pickle.load(handle2)

				with open(Input_directory + '/reporter_rev_avg_read__' + sample_name + '.txt', 'rb') as handle3:
					reporter_rev_avg_read_RNA = pickle.load(handle3)
					
				with open(Input_directory + '/genome_fwd_avg_read__' + sample_name + '.txt', 'rb') as handle4:
					genome_fwd_avg_read_RNA = pickle.load(handle4)

				with open(Input_directory + '/genome_rev_avg_read__' + sample_name + '.txt', 'rb') as handle5:
					genome_rev_avg_read_RNA = pickle.load(handle5)
					
				circuit_borders = (0, circuit_plasmid_length)
				reporter_borders = (0, reporter_plasmid_length)
				genome_borders = (0, genome_length)
				
				# ----------------------------------- #
				
				borders = circuit_borders
				
				profile_Forward_circuit = np.ones(borders[1]-borders[0])*1e-6
				for item in range(borders[0],borders[1]):
					for ss in range(len(circuit_fwd_avg_read_RNA[item])):
						if circuit_fwd_avg_read_RNA[item][ss] > 1:
							if item-borders[0]+circuit_fwd_avg_read_RNA[item][ss] < borders[1]-borders[0]:
								profile_Forward_circuit[item-borders[0]:item-borders[0]+circuit_fwd_avg_read_RNA[item][ss]] += np.ones(circuit_fwd_avg_read_RNA[item][ss])
					
				profile_Reverse_circuit = np.ones(borders[1]-borders[0])*1e-6
				for item in range(borders[0],borders[1]):
					for ss in range(len(circuit_rev_avg_read_RNA[item])):
						if circuit_rev_avg_read_RNA[item][ss] > 1:
							if item-borders[0]-circuit_rev_avg_read_RNA[item][ss]+1 >= 0:
								profile_Reverse_circuit[item-borders[0]-circuit_rev_avg_read_RNA[item][ss]+1:item-borders[0]+1] += np.ones(circuit_rev_avg_read_RNA[item][ss])

				Depth_of_fwd_read_circuit = sum(profile_Forward_circuit)
				Depth_of_rev_read_circuit = sum(profile_Reverse_circuit)
				
				# ----------------------------------- #

				borders = reporter_borders
			
				profile_Forward_reporter = np.ones(borders[1]-borders[0])*1e-6
				for item in range(borders[0],borders[1]):
					for ss in range(len(reporter_fwd_avg_read_RNA[item])):
						if reporter_fwd_avg_read_RNA[item][ss] > 1:
							if item-borders[0]+reporter_fwd_avg_read_RNA[item][ss] < borders[1]-borders[0]:
								profile_Forward_reporter[item-borders[0]:item-borders[0]+reporter_fwd_avg_read_RNA[item][ss]] += np.ones(reporter_fwd_avg_read_RNA[item][ss])
					
				profile_Reverse_reporter = np.ones(borders[1]-borders[0])*1e-6
				for item in range(borders[0],borders[1]):
					for ss in range(len(reporter_rev_avg_read_RNA[item])):
						if reporter_rev_avg_read_RNA[item][ss] > 1:
							if item-borders[0]-reporter_rev_avg_read_RNA[item][ss]+1 >= 0:
								profile_Reverse_reporter[item-borders[0]-reporter_rev_avg_read_RNA[item][ss]+1:item-borders[0]+1] += np.ones(reporter_rev_avg_read_RNA[item][ss])

				Depth_of_fwd_read_reporter = sum(profile_Forward_reporter)
				Depth_of_rev_read_reporter = sum(profile_Reverse_reporter)
				
				# ----------------------------------- #

				borders = genome_borders
			
				profile_Forward_genome = np.ones(borders[1]-borders[0])*1e-6
				for item in range(borders[0],borders[1]):
					for ss in range(len(genome_fwd_avg_read_RNA[item])):
						if genome_fwd_avg_read_RNA[item][ss] > 1:
							if item-borders[0]+genome_fwd_avg_read_RNA[item][ss] < borders[1]-borders[0]:
								profile_Forward_genome[item-borders[0]:item-borders[0]+genome_fwd_avg_read_RNA[item][ss]] += np.ones(genome_fwd_avg_read_RNA[item][ss])
					
				profile_Reverse_genome = np.ones(borders[1]-borders[0])*1e-6
				for item in range(borders[0],borders[1]):
					for ss in range(len(genome_rev_avg_read_RNA[item])):
						if genome_rev_avg_read_RNA[item][ss] > 1:
							if item-borders[0]-genome_rev_avg_read_RNA[item][ss]+1 >= 0:
								profile_Reverse_genome[item-borders[0]-genome_rev_avg_read_RNA[item][ss]+1:item-borders[0]+1] += np.ones(genome_rev_avg_read_RNA[item][ss])

				Depth_of_fwd_read_genome = sum(profile_Forward_genome)
				Depth_of_rev_read_genome = sum(profile_Reverse_genome)
				
				# ----------------------------------- #

				total_read = Depth_of_fwd_read_circuit + Depth_of_fwd_read_reporter + Depth_of_fwd_read_genome + Depth_of_rev_read_circuit + Depth_of_rev_read_reporter + Depth_of_rev_read_genome
				
				# ----------------------------------- #

				" deleting the profiles of regions within the tRNA and rRNA operons"

				profile_Forward_genome_tRNA_rRNA_operon_deleted = copy.deepcopy(profile_Forward_genome)
				profile_Reverse_genome_tRNA_rRNA_operon_deleted = copy.deepcopy(profile_Reverse_genome)

				for itemx in tRAN_rRNA_operon_range:
					profile_Forward_genome_tRNA_rRNA_operon_deleted[itemx[0]:itemx[1]] = np.ones(itemx[1]-itemx[0])*1e-6
					profile_Reverse_genome_tRNA_rRNA_operon_deleted[itemx[0]:itemx[1]] = np.ones(itemx[1]-itemx[0])*1e-6

				Depth_of_fwd_read_genome_tRNA_rRNA_operon_deleted = sum(profile_Forward_genome_tRNA_rRNA_operon_deleted)
				Depth_of_rev_read_genome_tRNA_rRNA_operon_deleted = sum(profile_Reverse_genome_tRNA_rRNA_operon_deleted)
					
				# ----------------------------------- #
				
				total_read_tRNA_rRNA_operon_deleted = Depth_of_fwd_read_circuit + Depth_of_fwd_read_reporter + Depth_of_fwd_read_genome_tRNA_rRNA_operon_deleted + Depth_of_rev_read_circuit + Depth_of_rev_read_reporter + Depth_of_rev_read_genome_tRNA_rRNA_operon_deleted
				
				# ----------------------------------- #

				file = open(self.output_RNA_seq_Profiles_filename + "/RNA_Seq_Profile__Circuit_Plasmid__" + sample_name + ".wig","w")
				
				file.write("Chromosome = " + circuit_plasmid_name + "\n")
				
				file.write("Normalized_transcription_profile_without_tRNA_rRNA_operons" + "\n")
				
				file.write("Position" + "\t" + "Forward_strand" + "\t" + "Reverse_strand" + "\n")
				
				for jj in range(len(profile_Forward_circuit)):
					file.write(str(jj+1) + "\t" + str((profile_Forward_circuit[jj]/total_read_tRNA_rRNA_operon_deleted)*1e9) + "\t" + str((profile_Reverse_circuit[jj]/total_read_tRNA_rRNA_operon_deleted)*1e9) + "\n")
				
				file.close()
				
				# ----------------------------------- #

				file = open(self.output_RNA_seq_Profiles_filename + "/RNA_Seq_Profile__Reporter_Plasmid__" + sample_name + ".wig","w")
				
				file.write("Chromosome = " + reporter_plasmid_name + "\n")
				
				file.write("Normalized_transcription_profile_without_tRNA_rRNA_operons" + "\n")
				
				file.write("Position" + "\t" + "Forward_strand" + "\t" + "Reverse_strand" + "\n")
				
				for jj in range(len(profile_Forward_reporter)):
					file.write(str(jj+1) + "\t" + str((profile_Forward_reporter[jj]/total_read_tRNA_rRNA_operon_deleted)*1e9) + "\t" + str((profile_Reverse_reporter[jj]/total_read_tRNA_rRNA_operon_deleted)*1e9) + "\n")
				
				file.close()
				
				# ----------------------------------- #

				file = open(self.output_RNA_seq_Profiles_filename + "/RNA_Seq_Profile__DH10B_Genome__" + sample_name + ".wig","w")
				
				file.write("Chromosome = " + genome_name + "\n")
				
				file.write("Normalized_transcription_profile_without_tRNA_rRNA_operons" + "\n")
				
				file.write("Position" + "\t" + "Forward_strand" + "\t" + "Reverse_strand" + "\n")
				
				for jj in range(len(profile_Forward_genome_tRNA_rRNA_operon_deleted)):
					file.write(str(jj+1) + "\t" + str((profile_Forward_genome_tRNA_rRNA_operon_deleted[jj]/total_read_tRNA_rRNA_operon_deleted)*1e9) + "\t" + str((profile_Reverse_genome_tRNA_rRNA_operon_deleted[jj]/total_read_tRNA_rRNA_operon_deleted)*1e9) + "\n")
				
				file.close()
				
				# ----------------------------------- #
				
				total_tRNA_fraction = (total_read - total_read_tRNA_rRNA_operon_deleted)/float(total_read)
				
				RNA_seq_total_read[sample_name] = total_read
				RNA_seq_total_read_tRNA_rRNA_operon_deleted[sample_name] = total_read_tRNA_rRNA_operon_deleted
				RNA_seq_total_tRNA_fraction[sample_name] = total_tRNA_fraction
				
				# ----------------------------------- #
				
				RNA_seq_Profile[sample_name] = {circuit_plasmid_name:{},reporter_plasmid_name:{},genome_name:{}}
				
				RNA_seq_Profile[sample_name][circuit_plasmid_name]['fwd'] = [(profile_Forward_circuit[jj]/total_read_tRNA_rRNA_operon_deleted)*1e9 for jj in range(len(profile_Forward_circuit)]
				RNA_seq_Profile[sample_name][circuit_plasmid_name]['rev'] = [(profile_Reverse_circuit[jj]/total_read_tRNA_rRNA_operon_deleted)*1e9 for jj in range(len(profile_Reverse_circuit)]
				
				RNA_seq_Profile[sample_name][reporter_plasmid_name]['fwd'] = [(profile_Forward_reporter[jj]/total_read_tRNA_rRNA_operon_deleted)*1e9 for jj in range(len(profile_Forward_reporter)]
				RNA_seq_Profile[sample_name][reporter_plasmid_name]['rev'] = [(profile_Reverse_reporter[jj]/total_read_tRNA_rRNA_operon_deleted)*1e9 for jj in range(len(profile_Reverse_reporter)]
				
				RNA_seq_Profile[sample_name][genome_name]['fwd'] = [(profile_Forward_genome_tRNA_rRNA_operon_deleted[jj]/total_read_tRNA_rRNA_operon_deleted)*1e9 for jj in range(len(profile_Forward_genome_tRNA_rRNA_operon_deleted)]
				RNA_seq_Profile[sample_name][genome_name]['rev'] = [(profile_Reverse_genome_tRNA_rRNA_operon_deleted[jj]/total_read_tRNA_rRNA_operon_deleted)*1e9 for jj in range(len(profile_Reverse_genome_tRNA_rRNA_operon_deleted)]
				
				
		return RNA_seq_Profile, RNA_seq_total_read, RNA_seq_total_read_tRNA_rRNA_operon_deleted, RNA_seq_total_tRNA_fraction
			
	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def generate_Ribo_Seq_profile(self):

		Ribo_seq_total_read = {}
		Ribo_seq_total_read_tRNA_rRNA_operon_deleted = {}
		Ribo_seq_total_tRNA_fraction = {}
		
		Ribo_seq_Profile = {}
		
		for sample_name in self.Coding_Genes_MAP_output_filename:
			
			if sample_name.split('_')[0] == 'RiboSeq':
				
				Input_directory = self.output_RiboSeq_reads_filename
			
				(circuit_plasmid_name, reporter_plasmid_name, genome_name) = self.reference_name[sample_name]
				(circuit_plasmid_length, reporter_plasmid_length, genome_length) = self.reference_length[sample_name]
				
				Master_wb_read = xlrd.open_workbook(self.tRNA_rRNA_operon_GFF_filename + '/NC_010473_1_operon.xls',"rb")
				Master_sh = Master_wb_read.sheet_by_name('NC_010473.1.operon')
				
				tRAN_rRNA_operon_range = []
				for i in range(43):
					left_bd = int(Master_sh.cell(2*i,11).value)
					right_bd = int(Master_sh.cell(2*i,12).value)
					tRAN_rRNA_operon_range.append((left_bd,right_bd))
				
				# ----------------------------------- #
				
				with open(Input_directory + '/circuit_fwd_avg_read__' + sample_name + '.txt', 'rb') as handle:
					circuit_fwd_avg_read_Ribo = pickle.load(handle)

				with open(Input_directory + '/circuit_rev_avg_read__' + sample_name + '.txt', 'rb') as handle1:
					circuit_rev_avg_read_Ribo = pickle.load(handle1)

				with open(Input_directory + '/reporter_fwd_avg_read__' + sample_name + '.txt', 'rb') as handle2:
					reporter_fwd_avg_read_Ribo = pickle.load(handle2)

				with open(Input_directory + '/reporter_rev_avg_read__' + sample_name + '.txt', 'rb') as handle3:
					reporter_rev_avg_read_Ribo = pickle.load(handle3)
					
				with open(Input_directory + '/genome_fwd_avg_read__' + sample_name + '.txt', 'rb') as handle4:
					genome_fwd_avg_read_Ribo = pickle.load(handle4)

				with open(Input_directory + '/genome_rev_avg_read__' + sample_name + '.txt', 'rb') as handle5:
					genome_rev_avg_read_Ribo = pickle.load(handle5)
					
				circuit_borders = (0, circuit_plasmid_length)
				reporter_borders = (0, reporter_plasmid_length)
				genome_borders = (0, genome_length)
				
				# ----------------------------------- #
				
				borders = circuit_borders
				
				profile_Forward_circuit = np.ones(borders[1]-borders[0])*1e-6
				for item in range(borders[0],borders[1]):
					for ss in range(len(circuit_fwd_avg_read_Ribo[item])):
						if circuit_fwd_avg_read_Ribo[item][ss] > 1:
							if item-borders[0]+circuit_fwd_avg_read_Ribo[item][ss] < borders[1]-borders[0] and circuit_fwd_avg_read_Ribo[item][ss] > 22:
								profile_Forward_circuit[item-borders[0]+11:item-borders[0]+circuit_fwd_avg_read_Ribo[item][ss]-11] += np.ones(circuit_fwd_avg_read_Ribo[item][ss]-22)/(circuit_fwd_avg_read_Ribo[item][ss]-22)
				
				profile_Reverse_circuit = np.ones(borders[1]-borders[0])*1e-6
				for item in range(borders[0],borders[1]):
					for ss in range(len(circuit_rev_avg_read_Ribo[item])):
						if circuit_rev_avg_read_Ribo[item][ss] > 1:
							if item-borders[0]-circuit_rev_avg_read_Ribo[item][ss] +1 >= 0 and circuit_rev_avg_read_Ribo[item][ss] > 22:
								profile_Reverse_circuit[item-borders[0]-circuit_rev_avg_read_Ribo[item][ss]+1+11:item-borders[0]+1-11] += np.ones(circuit_rev_avg_read_Ribo[item][ss]-22)/(circuit_rev_avg_read_Ribo[item][ss]-22)

				Depth_of_fwd_read_circuit = sum(profile_Forward_circuit)
				Depth_of_rev_read_circuit = sum(profile_Reverse_circuit)
				
				# ----------------------------------- #

				borders = reporter_borders
			
				profile_Forward_reporter = np.ones(borders[1]-borders[0])*1e-6
				for item in range(borders[0],borders[1]):
					for ss in range(len(reporter_fwd_avg_read_Ribo[item])):
						if reporter_fwd_avg_read_Ribo[item][ss] > 1:
							if item-borders[0]+reporter_fwd_avg_read_Ribo[item][ss] < borders[1]-borders[0] and reporter_fwd_avg_read_Ribo[item][ss] > 22:
								profile_Forward_reporter[item-borders[0]+11:item-borders[0]+reporter_fwd_avg_read_Ribo[item][ss]-11] += np.ones(reporter_fwd_avg_read_Ribo[item][ss]-22)/(reporter_fwd_avg_read_Ribo[item][ss]-22)
					
				profile_Reverse_reporter = np.ones(borders[1]-borders[0])*1e-6
				for item in range(borders[0],borders[1]):
					for ss in range(len(reporter_rev_avg_read_Ribo[item])):
						if reporter_rev_avg_read_Ribo[item][ss] > 1:
							if item-borders[0]-reporter_rev_avg_read_Ribo[item][ss]+1 >= 0 and reporter_rev_avg_read_Ribo[item][ss] > 22:
								profile_Reverse_reporter[item-borders[0]-reporter_rev_avg_read_Ribo[item][ss]+1+11:item-borders[0]+1-11] += np.ones(reporter_rev_avg_read_Ribo[item][ss]-22)/(reporter_rev_avg_read_Ribo[item][ss]-22)

				Depth_of_fwd_read_reporter = sum(profile_Forward_reporter)
				Depth_of_rev_read_reporter = sum(profile_Reverse_reporter)
				
				# ----------------------------------- #

				borders = genome_borders
			
				profile_Forward_genome = np.ones(borders[1]-borders[0])*1e-6
				for item in range(borders[0],borders[1]):
					for ss in range(len(genome_fwd_avg_read_Ribo[item])):
						if genome_fwd_avg_read_Ribo[item][ss] > 1:
							if item-borders[0]+genome_fwd_avg_read_Ribo[item][ss] < borders[1]-borders[0] and genome_fwd_avg_read_Ribo[item][ss] > 22:
								profile_Forward_genome[item-borders[0]+11:item-borders[0]+genome_fwd_avg_read_Ribo[item][ss]-11] += np.ones(genome_fwd_avg_read_Ribo[item][ss]-22)/(genome_fwd_avg_read_Ribo[item][ss]-22)
					
				profile_Reverse_genome = np.ones(borders[1]-borders[0])*1e-6
				for item in range(borders[0],borders[1]):
					for ss in range(len(genome_rev_avg_read_Ribo[item])):
						if genome_rev_avg_read_Ribo[item][ss] > 1:
							if item-borders[0]-genome_rev_avg_read_Ribo[item][ss]+1 >= 0 and genome_rev_avg_read_Ribo[item][ss] > 22:
								profile_Reverse_genome[item-borders[0]-genome_rev_avg_read_Ribo[item][ss]+1+11:item-borders[0]+1-11] += np.ones(genome_rev_avg_read_Ribo[item][ss]-22)/(genome_rev_avg_read_Ribo[item][ss]-22)

				Depth_of_fwd_read_genome = sum(profile_Forward_genome)
				Depth_of_rev_read_genome = sum(profile_Reverse_genome)
				
				# ----------------------------------- #

				total_read = Depth_of_fwd_read_circuit + Depth_of_fwd_read_reporter + Depth_of_fwd_read_genome + Depth_of_rev_read_circuit + Depth_of_rev_read_reporter + Depth_of_rev_read_genome
				
				# ----------------------------------- #

				" deleting the profiles of regions within the tRNA and rRNA operons"

				profile_Forward_genome_tRNA_rRNA_operon_deleted = copy.deepcopy(profile_Forward_genome)
				profile_Reverse_genome_tRNA_rRNA_operon_deleted = copy.deepcopy(profile_Reverse_genome)

				for itemx in tRAN_rRNA_operon_range:
					profile_Forward_genome_tRNA_rRNA_operon_deleted[itemx[0]:itemx[1]] = np.ones(itemx[1]-itemx[0])*1e-6
					profile_Reverse_genome_tRNA_rRNA_operon_deleted[itemx[0]:itemx[1]] = np.ones(itemx[1]-itemx[0])*1e-6

				Depth_of_fwd_read_genome_tRNA_rRNA_operon_deleted = sum(profile_Forward_genome_tRNA_rRNA_operon_deleted)
				Depth_of_rev_read_genome_tRNA_rRNA_operon_deleted = sum(profile_Reverse_genome_tRNA_rRNA_operon_deleted)
					
				# ----------------------------------- #
				
				total_read_tRNA_rRNA_operon_deleted = Depth_of_fwd_read_circuit + Depth_of_fwd_read_reporter + Depth_of_fwd_read_genome_tRNA_rRNA_operon_deleted + Depth_of_rev_read_circuit + Depth_of_rev_read_reporter + Depth_of_rev_read_genome_tRNA_rRNA_operon_deleted
				
				# ----------------------------------- #

				file = open(self.output_Ribo_seq_Profiles_filename + "/Ribosome_Occupancy_Profile__Circuit_Plasmid__" + sample_name + ".wig","w")
				
				file.write("Chromosome = " + circuit_plasmid_name + "\n")
				
				file.write("Ribosome_occupancy" + "\n")
				
				file.write("Position" + "\t" + "Forward_strand" + "\t" + "Reverse_strand" + "\n")
				
				for jj in range(len(profile_Forward_circuit)):
					#file.write(str(jj+1) + "\t" + str((profile_Forward_circuit[jj]/total_read_tRNA_rRNA_operon_deleted)*self.total_active_ribosomes) + "\t" + str((profile_Reverse_circuit[jj]/total_read_tRNA_rRNA_operon_deleted)*self.total_active_ribosomes) + "\n")
					file.write(str(jj+1) + "\t" + str((profile_Forward_circuit[jj]/total_read)*self.total_active_ribosomes) + "\t" + str((profile_Reverse_circuit[jj]/total_read)*self.total_active_ribosomes) + "\n")

				file.close()
				
				# ----------------------------------- #

				file = open(self.output_Ribo_seq_Profiles_filename + "/Ribosome_Occupancy_Profile__Reporter_Plasmid__" + sample_name + ".wig","w")
				
				file.write("Chromosome = " + reporter_plasmid_name + "\n")
				
				file.write("Ribosome_occupancy" + "\n")
				
				file.write("Position" + "\t" + "Forward_strand" + "\t" + "Reverse_strand" + "\n")
				
				for jj in range(len(profile_Forward_reporter)):
					#file.write(str(jj+1) + "\t" + str((profile_Forward_reporter[jj]/total_read_tRNA_rRNA_operon_deleted)*self.total_active_ribosomes) + "\t" + str((profile_Reverse_reporter[jj]/total_read_tRNA_rRNA_operon_deleted)*self.total_active_ribosomes) + "\n")
					file.write(str(jj+1) + "\t" + str((profile_Forward_reporter[jj]/total_read)*self.total_active_ribosomes) + "\t" + str((profile_Reverse_reporter[jj]/total_read)*self.total_active_ribosomes) + "\n")
				
				file.close()
				
				# ----------------------------------- #

				file = open(self.output_Ribo_seq_Profiles_filename + "/Ribosome_Occupancy_Profile__DH10B_Genome__" + sample_name + ".wig","w")
				
				file.write("Chromosome = " + genome_name + "\n")
				
				file.write("Ribosome_occupancy" + "\n")
				
				file.write("Position" + "\t" + "Forward_strand" + "\t" + "Reverse_strand" + "\n")
				
				for jj in range(len(profile_Forward_genome)):
					#file.write(str(jj+1) + "\t" + str((profile_Forward_genome_tRNA_rRNA_operon_deleted[jj]/total_read_tRNA_rRNA_operon_deleted)*self.total_active_ribosomes) + "\t" + str((profile_Reverse_genome_tRNA_rRNA_operon_deleted[jj]/total_read_tRNA_rRNA_operon_deleted)*self.total_active_ribosomes) + "\n")
					file.write(str(jj+1) + "\t" + str((profile_Forward_genome[jj]/total_read)*self.total_active_ribosomes) + "\t" + str((profile_Reverse_genome[jj]/total_read)*self.total_active_ribosomes) + "\n")
				
				file.close()
				
				# ----------------------------------- #
				
				total_tRNA_fraction = (total_read - total_read_tRNA_rRNA_operon_deleted)/float(total_read)
		
				Ribo_seq_total_read[sample_name] = total_read
				Ribo_seq_total_read_tRNA_rRNA_operon_deleted[sample_name] = total_read_tRNA_rRNA_operon_deleted
				Ribo_seq_total_tRNA_fraction[sample_name] = total_tRNA_fraction
		
				# ----------------------------------- #
				
				Ribo_seq_Profile[sample_name] = {circuit_plasmid_name:{},reporter_plasmid_name:{},genome_name:{}}
				
				#Ribo_seq_Profile[sample_name][circuit_plasmid_name]['fwd'] = [(profile_Forward_circuit[jj]/total_read_tRNA_rRNA_operon_deleted)*self.total_active_ribosomes for jj in range(len(profile_Forward_circuit)]
				#Ribo_seq_Profile[sample_name][circuit_plasmid_name]['rev'] = [(profile_Reverse_circuit[jj]/total_read_tRNA_rRNA_operon_deleted)*self.total_active_ribosomes for jj in range(len(profile_Reverse_circuit)]
				
				#Ribo_seq_Profile[sample_name][reporter_plasmid_name]['fwd'] = [(profile_Forward_reporter[jj]/total_read_tRNA_rRNA_operon_deleted)*self.total_active_ribosomes for jj in range(len(profile_Forward_reporter)]
				#Ribo_seq_Profile[sample_name][reporter_plasmid_name]['rev'] = [(profile_Reverse_reporter[jj]/total_read_tRNA_rRNA_operon_deleted)*self.total_active_ribosomes for jj in range(len(profile_Reverse_reporter)]
				
				#Ribo_seq_Profile[sample_name][genome_name]['fwd'] = [(profile_Forward_genome_tRNA_rRNA_operon_deleted[jj]/total_read_tRNA_rRNA_operon_deleted)*self.total_active_ribosomes for jj in range(len(profile_Forward_genome_tRNA_rRNA_operon_deleted)]
				#Ribo_seq_Profile[sample_name][genome_name]['rev'] = [(profile_Reverse_genome_tRNA_rRNA_operon_deleted[jj]/total_read_tRNA_rRNA_operon_deleted)*self.total_active_ribosomes for jj in range(len(profile_Reverse_genome_tRNA_rRNA_operon_deleted)]
				
				
				Ribo_seq_Profile[sample_name][circuit_plasmid_name]['fwd'] = [(profile_Forward_circuit[jj]/total_read)*self.total_active_ribosomes for jj in range(len(profile_Forward_circuit)]
				Ribo_seq_Profile[sample_name][circuit_plasmid_name]['rev'] = [(profile_Reverse_circuit[jj]/total_read)*self.total_active_ribosomes for jj in range(len(profile_Reverse_circuit)]
				
				Ribo_seq_Profile[sample_name][reporter_plasmid_name]['fwd'] = [(profile_Forward_reporter[jj]/total_read)*self.total_active_ribosomes for jj in range(len(profile_Forward_reporter)]
				Ribo_seq_Profile[sample_name][reporter_plasmid_name]['rev'] = [(profile_Reverse_reporter[jj]/total_read)*self.total_active_ribosomes for jj in range(len(profile_Reverse_reporter)]
				
				Ribo_seq_Profile[sample_name][genome_name]['fwd'] = [(profile_Forward_genome[jj]/total_read)*self.total_active_ribosomes for jj in range(len(profile_Forward_genome)]
				Ribo_seq_Profile[sample_name][genome_name]['rev'] = [(profile_Reverse_genome[jj]/total_read)*self.total_active_ribosomes for jj in range(len(profile_Reverse_genome)]
		
		
		return Ribo_seq_Profile, Ribo_seq_total_read, Ribo_seq_total_read_tRNA_rRNA_operon_deleted, Ribo_seq_total_tRNA_fraction

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def generate_mapped_fragment_size_distribution(self):

		fragment_size_distribution = {}
		
		for sample_name in self.Coding_Genes_MAP_output_filename:
			
			
			if sample_name.split('_')[0] == 'RNASeq':
				Input_directory = self.output_RNASeq_reads_filename
			elif sample_name.split('_')[0] == 'RiboSeq':
				Input_directory = self.output_RiboSeq_reads_filename
			
			
			with open(Input_directory + '/circuit_fwd_avg_read__' + sample_name + '.txt', 'rb') as handle:
				circuit_fwd_avg_read = pickle.load(handle)

			with open(Input_directory + '/circuit_rev_avg_read__' + sample_name + '.txt', 'rb') as handle1:
				circuit_rev_avg_read = pickle.load(handle1)

			with open(Input_directory + '/reporter_fwd_avg_read__' + sample_name + '.txt', 'rb') as handle2:
				reporter_fwd_avg_read = pickle.load(handle2)

			with open(Input_directory + '/reporter_rev_avg_read__' + sample_name + '.txt', 'rb') as handle3:
				reporter_rev_avg_read = pickle.load(handle3)
				
			with open(Input_directory + '/genome_fwd_avg_read__' + sample_name + '.txt', 'rb') as handle4:
				genome_fwd_avg_read = pickle.load(handle4)

			with open(Input_directory + '/genome_rev_avg_read__' + sample_name + '.txt', 'rb') as handle5:
				genome_rev_avg_read = pickle.load(handle5)
				
			# ----------------------------------- #
			
			all_frag_plasmid = []
			
			for ppp in circuit_fwd_avg_read.values():
				all_frag_plasmid += ppp[1:len(ppp)]
					
			for aaa in circuit_rev_avg_read.values():
				all_frag_plasmid += aaa[1:len(aaa)]

			for ppp in reporter_fwd_avg_read.values():
				all_frag_plasmid += ppp[1:len(ppp)]
					
			for aaa in reporter_rev_avg_read.values():
				all_frag_plasmid += aaa[1:len(aaa)]

			for ppp in genome_fwd_avg_read.values():
				all_frag_plasmid += ppp[1:len(ppp)]
					
			for aaa in genome_rev_avg_read.values():
				all_frag_plasmid += aaa[1:len(aaa)]

			flag = True
			while flag:
				if all_frag_plasmid.count(1) > 0:
					all_frag_plasmid.remove(1)
				else:
					flag = False
			
			fragdist = np.array(range(self.frag_dist_bin[0],self.frag_dist_bin[1]))
			frag_dist_probability = np.histogram(all_frag_plasmid, bins=fragdist, density=True)
			
			fragment_size_distribution[sample_name] = (frag_dist_probability[1][1:], frag_dist_probability[0])
			
		return fragment_size_distribution

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def identify_promoter_TSS(self, RNA_seq_Profile, sample_name, Chromosome, Borders, alfa):
		
		profile = RNA_seq_Profile[sample_name][Chromosome]['fwd'][Borders[0]:Borders[1]]
		RNAP_Flux_Profile = [alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile]
		min_value = min([x for x in RNAP_Flux_Profile if x > 7.2511e-10])
		diff_log = np.diff(np.log10([max(min_value, x) for x in RNAP_Flux_Profile]))
		diff_log_pos = [max(0, x) for x in diff_log]
		
		profile_rev = RNA_seq_Profile[sample_name][Chromosome]['rev'][Borders[0]:Borders[1]]
		RNAP_Flux_Profile_rev = [alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile_rev]
		min_value_rev = min([x for x in RNAP_Flux_Profile_rev if x > 7.2511e-10])
		diff_log_rev = np.diff(np.log10([max(min_value_rev, x) for x in RNAP_Flux_Profile_rev]))
		diff_log_neg_rev = [min(0, x) for x in diff_log_rev]
		
		TSS_fwd = {sample_name:{Chromosome:[]}}
		TSS_rev = {sample_name:{Chromosome:[]}}
		
		for pos in range(len(diff_log_pos)):
			if diff_log_pos[pos] > 0.7:
				TSS_fwd[sample_name][Chromosome].append(Borders[0]+pos)
				
		for pos in range(len(diff_log_neg_rev)):
			if diff_log_neg_rev[pos] < -0.7:
				TSS_rev[sample_name][Chromosome].append(Borders[0]+pos)

		return diff_log_pos, diff_log_neg_rev, TSS_fwd, TSS_rev

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def calculate_promoter_activity(self, RNA_seq_Profile, sample_name, Chromosome, Strand, TSS, X2, X1, alfa):

		Borders = [TSS-X2,
			TSS-X1,
			TSS+X1,
			TSS+X2]

		profile_upstream = RNA_seq_Profile[sample_name][Chromosome][Strand][Borders[0]:Borders[1]]
		UP = np.mean([alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile_upstream])

		profile_downstream = RNA_seq_Profile[sample_name][Chromosome][Strand][Borders[2]:Borders[3]]
		DOWN = np.mean([alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile_downstream])
				
		if Strand == 'fwd':
			promoter_activity = max(2e-7, DOWN - UP)
		elif Strand == 'rev':
			promoter_activity = max(2e-7, UP - DOWN)
			
		return promoter_activity

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def calculate_Cello_promoter_activity(self, Promoter_list, Promoter_list_2, alfa, alfa2):

		Cello_flux = {item:[] for item in Promoter_list}
		Cello_flux2 = {item:[] for item in Promoter_list_2}

		for uu in range(len(self.inducer)):
			print "-------------state-------------", uu
			
			flux = {}

			flux['pBAD1']    = (self.par['pBAD1'][self.inducer[uu+1]['pBAD1']])
			flux['pTet1']    = (self.par['pTet1'][self.inducer[uu+1]['pTet1']])
			flux['SrpR']     = (self.par['pBAD1'][self.inducer[uu+1]['pBAD1']] + self.par['pTet1'][self.inducer[uu+1]['pTet1']])

			flux['pBAD2']    = (self.par['pBAD2'][self.inducer[uu+1]['pBAD2']])
			flux['AmtR']     = (self.par['pBAD2'][self.inducer[uu+1]['pBAD2']])

			flux['pTet2']    = (self.par['pTet2'][self.inducer[uu+1]['pTet2']])
			flux['AmeR']     = (self.par['pTet2'][self.inducer[uu+1]['pTet2']])

			flux['pAmtR']    = (self.par['pAmtR']['min']+(self.par['pAmtR']['max']-self.par['pAmtR']['min'])*((self.par['pAmtR']['K']**self.par['pAmtR']['n'])/((self.par['pAmtR']['K']**self.par['pAmtR']['n'])+(flux['AmtR']**self.par['pAmtR']['n']))))
			flux['pAmeR']    = (self.par['pAmeR']['min']+(self.par['pAmeR']['max']-self.par['pAmeR']['min'])*((self.par['pAmeR']['K']**self.par['pAmeR']['n'])/((self.par['pAmeR']['K']**self.par['pAmeR']['n'])+(flux['AmeR']**self.par['pAmeR']['n']))))
			flux['BetI']     = (self.par['pAmtR']['min']+(self.par['pAmtR']['max']-self.par['pAmtR']['min'])*((self.par['pAmtR']['K']**self.par['pAmtR']['n'])/((self.par['pAmtR']['K']**self.par['pAmtR']['n'])+(flux['AmtR']**self.par['pAmtR']['n']))) + self.par['pAmeR']['min']+(self.par['pAmeR']['max']-self.par['pAmeR']['min'])*((self.par['pAmeR']['K']**self.par['pAmeR']['n'])/((self.par['pAmeR']['K']**self.par['pAmeR']['n'])+(flux['AmeR']**self.par['pAmeR']['n']))))

			flux['pSrpR']    = (self.par['pSrpR']['min']+(self.par['pSrpR']['max']-self.par['pSrpR']['min'])*((self.par['pSrpR']['K']**self.par['pSrpR']['n'])/((self.par['pSrpR']['K']**self.par['pSrpR']['n'])+(flux['SrpR']**self.par['pSrpR']['n']))))
			flux['pBetI']    = (self.par['pBetI']['min']+(self.par['pBetI']['max']-self.par['pBetI']['min'])*((self.par['pBetI']['K']**self.par['pBetI']['n'])/((self.par['pBetI']['K']**self.par['pBetI']['n'])+(flux['BetI']**self.par['pBetI']['n']))))
			flux['PhlF']     = (self.par['pSrpR']['min']+(self.par['pSrpR']['max']-self.par['pSrpR']['min'])*((self.par['pSrpR']['K']**self.par['pSrpR']['n'])/((self.par['pSrpR']['K']**self.par['pSrpR']['n'])+(flux['SrpR']**self.par['pSrpR']['n']))) + self.par['pBetI']['min']+(self.par['pBetI']['max']-self.par['pBetI']['min'])*((self.par['pBetI']['K']**self.par['pBetI']['n'])/((self.par['pBetI']['K']**self.par['pBetI']['n'])+(flux['BetI']**self.par['pBetI']['n']))))

			flux['pTac']     = (self.par['pTac'][self.inducer[uu+1]['pTac']])
			flux['HlyIIR']   = (self.par['pTac'][self.inducer[uu+1]['pTac']])

			flux['pPhlF']    = (self.par['pPhlF']['min']+(self.par['pPhlF']['max']-self.par['pPhlF']['min'])*((self.par['pPhlF']['K']**self.par['pPhlF']['n'])/((self.par['pPhlF']['K']**self.par['pPhlF']['n'])+(flux['PhlF']**self.par['pPhlF']['n']))))
			flux['pHlyIIR']  = (self.par['pHlyIIR']['min']+(self.par['pHlyIIR']['max']-self.par['pHlyIIR']['min'])*((self.par['pHlyIIR']['K']**self.par['pHlyIIR']['n'])/((self.par['pHlyIIR']['K']**self.par['pHlyIIR']['n'])+(flux['HlyIIR']**self.par['pHlyIIR']['n']))))
			flux['BM3R1']    = (self.par['pPhlF']['min']+(self.par['pPhlF']['max']-self.par['pPhlF']['min'])*((self.par['pPhlF']['K']**self.par['pPhlF']['n'])/((self.par['pPhlF']['K']**self.par['pPhlF']['n'])+(flux['PhlF']**self.par['pPhlF']['n']))) + self.par['pHlyIIR']['min']+(self.par['pHlyIIR']['max']-self.par['pHlyIIR']['min'])*((self.par['pHlyIIR']['K']**self.par['pHlyIIR']['n'])/((self.par['pHlyIIR']['K']**self.par['pHlyIIR']['n'])+(flux['HlyIIR']**self.par['pHlyIIR']['n']))))

			flux['pBM3R1']  = (self.par['pBM3R1']['min']+(self.par['pBM3R1']['max']-self.par['pBM3R1']['min'])*((self.par['pBM3R1']['K']**self.par['pBM3R1']['n'])/((self.par['pBM3R1']['K']**self.par['pBM3R1']['n'])+(flux['BM3R1']**self.par['pBM3R1']['n']))))

			for item in Cello_flux:
				Cello_flux[item].append(alfa * self.RNAP_p15A * flux[item])		# unit RNAP/s
				
			for item in Cello_flux2:
				Cello_flux2[item].append(alfa2 * self.RNAP_p15A * flux[item])		# unit RNAP/s
			
		return Cello_flux, Cello_flux2

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def calculate_ribozyme_cleavege_efficiency(self, RNA_seq_Reads, sample_name, Chromosome, Strand, ribozyme_cut_position, offset):
		
		cut_left_profile = np.ones(2*offset)*1e-6
		uncut_profile = np.ones(2*offset)*1e-6
		cut_right_profile = np.ones(2*offset)*1e-6
		
		for pos in range(ribozyme_cut_position-offset,ribozyme_cut_position+offset):
			
			if pos < ribozyme_cut_position:
				
				for ss in range(len(RNA_seq_Reads[sample_name][Chromosome][Strand][pos])):
					
					read_length = RNA_seq_Reads[sample_name][Chromosome][Strand][pos][ss]
					
					if read_length > 1 and pos-(ribozyme_cut_position-offset)+read_length < 2*offset:
						
						if pos+read_length == ribozyme_cut_position-1:
							
							cut_left_profile[pos-(ribozyme_cut_position-offset):pos-(ribozyme_cut_position-offset)+read_length] += np.ones(read_length)
							
						elif pos+read_length > ribozyme_cut_position-1:
							
							uncut_profile[pos-(ribozyme_cut_position-offset):pos-(ribozyme_cut_position-offset)+read_length] += np.ones(read_length)
							
			elif pos == ribozyme_cut_position:
				
				for ss in range(len(RNA_seq_Reads[sample_name][Chromosome][Strand][pos])):
					
					read_length = RNA_seq_Reads[sample_name][Chromosome][Strand][pos][ss]
					
					if read_length > 1 and pos-(ribozyme_cut_position-offset)+read_length < 2*offset:
						
						cut_right_profile[pos-(ribozyme_cut_position-offset):pos-(ribozyme_cut_position-offset)+read_length] += np.ones(read_length)

		ribozyme_cleavege_efficiency = 1.0 -  float(sum(uncut_profile[offset:]))/(sum(cut_left_profile[offset:]) + sum(uncut_profile[offset:]) + sum(cut_right_profile[offset:]))
		
		return ribozyme_cleavege_efficiency

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def identify_terminator_TTS(self, RNA_seq_Profile, sample_name, Chromosome, Borders, alfa, window_size):


		profile = RNA_seq_Profile[sample_name][Chromosome]['fwd'][Borders[0]:Borders[1]]
		RNAP_Flux_Profile = [alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile]
		min_value = min([x for x in RNAP_Flux_Profile if x > 7.2511e-10])
		window_averaged_RNAP_Flux_Profile = [1 for kkk in range(window_size)] + [float(max(min_value,np.mean(RNAP_Flux_Profile[jjj-window_size:jjj]))) / max(min_value,np.mean(RNAP_Flux_Profile[jjj:jjj+window_size])) for jjj in range(window_size,len(RNAP_Flux_Profile)-window_size)] + [1 for kkk in range(window_size)]
		diff_log = np.log10(window_averaged_RNAP_Flux_Profile)
		diff_log_pos = [max(0, x) for x in diff_log]
		second_diff_log = np.diff([x if x > 0.7 else 0.0 for x in diff_log])
		
		TTS_fwd = {sample_name:{Chromosome:[]}}
		for wdw in range(1,len(second_diff_log)):
			if second_diff_log[wdw-1] > 0 and second_diff_log[wdw] <= 0:
				TTS_fwd[sample_name][Chromosome].append(Borders[0]+wdw)
		
		profile_rev = RNA_seq_Profile[sample_name][Chromosome]['rev'][Borders[0]:Borders[1]]
		RNAP_Flux_Profile_rev = [alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile_rev]
		min_value_rev = min([x for x in RNAP_Flux_Profile_rev if x > 7.2511e-10])
		window_averaged_RNAP_Flux_Profile_rev = [1 for kkk in range(window_size)] + [float(max(min_value_rev,np.mean(RNAP_Flux_Profile_rev[jjj-window_size:jjj]))) / max(min_value_rev,np.mean(RNAP_Flux_Profile_rev[jjj:jjj+window_size])) for jjj in range(window_size,len(RNAP_Flux_Profile_rev)-window_size)] + [1 for kkk in range(window_size)]
		diff_log_rev = np.log10(window_averaged_RNAP_Flux_Profile_rev)
		diff_log_neg_rev = [min(0, x) for x in diff_log_rev]
		second_diff_log_rev = np.diff([x if x < -0.7 else 0.0 for x in diff_log_rev])
		
		TTS_rev = {sample_name:{Chromosome:[]}}
		for wdw in range(1,len(second_diff_log_rev)):
			if second_diff_log_rev[wdw-1] < 0 and second_diff_log_rev[wdw] >= 0:
				TTS_rev[sample_name][Chromosome].append(Borders[0]+wdw)
		
		TTS_window_fwd = {sample_name:{Chromosome:[]}}
		TTS_window_rev = {sample_name:{Chromosome:[]}}
		
		for pos in range(len(diff_log_pos)):
			if diff_log_pos[pos] > 0.7:
				TTS_window_fwd[sample_name][Chromosome].append(Borders[0]+pos)
				
		for pos in range(len(diff_log_neg_rev)):
			if diff_log_neg_rev[pos] < -0.7:
				TTS_window_rev[sample_name][Chromosome].append(Borders[0]+pos)

		return diff_log_pos, diff_log_neg_rev, TTS_window_fwd, TTS_window_rev, second_diff_log, second_diff_log_rev, TTS_fwd, TTS_rev

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def calculate_terminator_strength(self, RNA_seq_Profile, sample_name, Chromosome, Strand, TTS, X2, X1, alfa):

		Borders = [TTS-X2,
			TTS-X1,
			TTS+X1,
			TTS+X2]

		profile_upstream = RNA_seq_Profile[sample_name][Chromosome][Strand][Borders[0]:Borders[1]]
		RNAP_Flux_Profile_up = [alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile_upstream]
		UP = np.mean(RNAP_Flux_Profile_up)

		profile_downstream = RNA_seq_Profile[sample_name][Chromosome][Strand][Borders[2]:Borders[3]]
		RNAP_Flux_Profile_down = [alfa * self.RNAP_p15A * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile_downstream]
		DOWN = np.mean(RNAP_Flux_Profile_down)
				
		if Strand == 'fwd':
			terminator_strength = float(UP)/DOWN
		elif Strand == 'rev':
			terminator_strength = float(DOWN)/UP
			
		return terminator_strength

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def calculate_FPKM(self, GFF_filename, RNA_seq_Profile, sample_name, Chromosome, gene):

		[circuit_plasmid, reporter_plasmid, genome] = self.reference_name[sample_name]
		
		plasmids_order = [circuit_plasmid, reporter_plasmid, genome]
		Gene_dict = self.read_gff([GFF_filename],plasmids_order)
		
		start = Gene_dict[Chromosome][gene]['start']
		end = Gene_dict[Chromosome][gene]['end']
		if Gene_dict[Chromosome][gene]['strand'] == '+':
			FPKM = np.mean(RNA_seq_Profile[sample_name][Chromosome]['fwd'][start:end])
		elif Gene_dict[Chromosome][gene]['strand'] == '-':
			FPKM = np.mean(RNA_seq_Profile[sample_name][Chromosome]['rev'][start:end])
		
		return FPKM
	
	# ---------------------------------------------------------------------------------------------------------------------------------------- #
	
	def generate_FPKM_table(self, GFF_filename, RNA_seq_Profile):
		
		FPKM_dict = {}
		
		for sample_name in self.reference_name:
			
			if sample_name.split('_')[0] == 'RNASeq':
			
				[circuit_plasmid, reporter_plasmid, genome] = self.reference_name[sample_name]
				
				plasmids_order = [circuit_plasmid, reporter_plasmid, genome]
				Gene_dict = self.read_gff([GFF_filename],plasmids_order)

				FPKM_dict[sample_name] = {circuit_plasmid:{},reporter_plasmid:{},genome:{}}
				
				for item in Gene_dict[circuit_plasmid]:
					if Gene_dict[circuit_plasmid][item]['feature'] in ['gene','resistance_maker']:
						start = Gene_dict[circuit_plasmid][item]['start']
						end = Gene_dict[circuit_plasmid][item]['end']
						if Gene_dict[circuit_plasmid][item]['strand'] == '+':
							FPKM_dict[sample_name][circuit_plasmid][item] = np.mean(RNA_seq_Profile[sample_name][circuit_plasmid]['fwd'][start:end])
						elif Gene_dict[circuit_plasmid][item]['strand'] == '-':
							FPKM_dict[sample_name][circuit_plasmid][item] = np.mean(RNA_seq_Profile[sample_name][circuit_plasmid]['rev'][start:end])
				
				for item in Gene_dict[reporter_plasmid]:
					if Gene_dict[reporter_plasmid][item]['feature'] in ['gene','resistance_maker']:
						start = Gene_dict[reporter_plasmid][item]['start']
						end = Gene_dict[reporter_plasmid][item]['end']
						if Gene_dict[reporter_plasmid][item]['strand'] == '+':
							FPKM_dict[sample_name][reporter_plasmid][item] = np.mean(RNA_seq_Profile[sample_name][reporter_plasmid]['fwd'][start:end])
						elif Gene_dict[reporter_plasmid][item]['strand'] == '-':
							FPKM_dict[sample_name][reporter_plasmid][item] = np.mean(RNA_seq_Profile[sample_name][reporter_plasmid]['rev'][start:end])
					
				for item in Gene_dict[genome]:
					start = Gene_dict[genome][item]['start']
					end = Gene_dict[genome][item]['end']
					if Gene_dict[genome][item]['strand'] == '+':
						FPKM_dict[sample_name][genome][item] = np.mean(RNA_seq_Profile[sample_name][genome]['fwd'][start:end])
					elif Gene_dict[genome][item]['strand'] == '-':
						FPKM_dict[sample_name][genome][item] = np.mean(RNA_seq_Profile[sample_name][genome]['rev'][start:end])

		# ----------------------------------- #

		workbook = xlwt.Workbook()

		worksheet = workbook.add_sheet('FPKM')

		worksheet.write(0,0,"Reference")
		worksheet.write(0,1,"Gene")
		
		[circuit_plasmid, reporter_plasmid, genome] = self.reference_name['RNASeq_State_1']
		
		row2 = 1
		for chromos in self.reference_name['RNASeq_State_1']:
			for gene2 in sorted(Gene_dict[chromos]):
				if Gene_dict[chromos][gene2]['feature'] in ['gene','resistance_maker']:
					if chromos == circuit_plasmid:
						worksheet.write(row2,0,'Circuit_Plasmid')
					elif chromos == reporter_plasmid:
						worksheet.write(row2,0,'Reporter_Plasmid')
					elif chromos == genome:
						worksheet.write(row2,0,'NC_010473.1')
					worksheet.write(row2,1,gene2.split()[0])
					row2 += 1
		
		kk = 2
		for sample_name in sorted(FPKM_dict):
			worksheet.write(0,kk,sample_name)
			row2 = 1
			for chromos in self.reference_name['RNASeq_State_1']:
				for gene2 in sorted(Gene_dict[chromos]):
					if Gene_dict[chromos][gene2]['feature'] in ['gene','resistance_maker']:
						if sample_name == 'RNASeq_Control':
							[circuit_plasmid2, reporter_plasmid2, genome2] = self.reference_name[sample_name]
							
							if chromos == circuit_plasmid:
								chromos2 = circuit_plasmid2
							elif chromos == reporter_plasmid:
								chromos2 = reporter_plasmid2
							else:
								chromos2 = chromos
							
							if gene2 in FPKM_dict[sample_name][chromos2]:
								worksheet.write(row2,kk,FPKM_dict[sample_name][chromos2][gene2])
							else:
								worksheet.write(row2,kk,0)
							row2 += 1
						else:
							worksheet.write(row2,kk,FPKM_dict[sample_name][chromos][gene2])
							row2 += 1
			
			kk += 1
			
		workbook.save(self.output_FPKM_table_filename+'/FPKM.xls')

		return FPKM_dict

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def calculate_mRNA_level(self, GFF_filename, RNA_seq_Profile, sample_name, Chromosome, gene, alfa):
		
		[circuit_plasmid, reporter_plasmid, genome] = self.reference_name[sample_name]
		
		plasmids_order = [circuit_plasmid, reporter_plasmid, genome]
		Gene_dict = self.read_gff([GFF_filename],plasmids_order)
		
		start = Gene_dict[Chromosome][gene]['start']
		end = Gene_dict[Chromosome][gene]['end']
		
		if Gene_dict[Chromosome][gene]['strand'] == '+':
			
			if gene == 'AmtR':
				Borders = (end-150-20,end-150-10)
			elif gene == 'SrpR':
				Borders = (end-20,end-10)
			else:
				Borders = (end-10,end-0)
			
			profile = RNA_seq_Profile[sample_name][Chromosome]['fwd'][Borders[0]:Borders[1]]
			RNAP_Flux = np.mean([alfa * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])
			
		elif Gene_dict[Chromosome][gene]['strand'] == '-':
			
			Borders = (start,start+10)
			
			profile = RNA_seq_Profile[sample_name][Chromosome]['rev'][Borders[0]:Borders[1]]
			RNAP_Flux = np.mean([alfa * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])
		
		mRNA = float(RNAP_Flux) / self.degradation

		return mRNA

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def calculate_translation_efficiency(self, GFF_filename, RNA_seq_Profile, sample_name, Chromosome, gene, alfa, RD_Data):
		
		[circuit_plasmid, reporter_plasmid, genome] = self.reference_name[sample_name]
		
		plasmids_order = [circuit_plasmid, reporter_plasmid, genome]
		Gene_dict = self.read_gff([GFF_filename],plasmids_order)
		
		start = Gene_dict[Chromosome][gene]['start']
		end = Gene_dict[Chromosome][gene]['end']
		
		if Gene_dict[Chromosome][gene]['strand'] == '+':
			
			if gene == 'AmtR':
				Borders = (end-150-20,end-150-10)
			elif gene == 'SrpR':
				Borders = (end-20,end-10)
			else:
				Borders = (end-10,end-0)
			
			profile = RNA_seq_Profile[sample_name][Chromosome]['fwd'][Borders[0]:Borders[1]]
			RNAP_Flux = np.mean([alfa * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])
			
		elif Gene_dict[Chromosome][gene]['strand'] == '-':
			
			Borders = (start,start+10)
			
			profile = RNA_seq_Profile[sample_name][Chromosome]['rev'][Borders[0]:Borders[1]]
			RNAP_Flux = np.mean([alfa * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])
		
		mRNA = float(RNAP_Flux) / self.degradation
		
		sample_name_Ribo = 'RiboSeq' + sample_name.split('RNASeq')[1]
		
		TE = float(RD_Data[sample_name_Ribo][Chromosome][gene] * self.omega)/mRNA
		
		return TE

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def generate_RD_table(self, GFF_filename, Ribo_seq_Profile, aa_ends_excluded):

		RD_dict = {}
		
		for sample_name in self.reference_name:

			if sample_name.split('_')[0] == 'RiboSeq':

				[circuit_plasmid, reporter_plasmid, genome] = self.reference_name[sample_name]
				
				plasmids_order = [circuit_plasmid, reporter_plasmid, genome]
				Gene_dict = self.read_gff([GFF_filename],plasmids_order)

				RD_dict[sample_name] = {circuit_plasmid:{},reporter_plasmid:{},genome:{}}
				
				# ----------------------------------- #
				
				circuit_gene_dict_P = []
				circuit_gene_dict_M = []
				for item in Gene_dict[circuit_plasmid]:
					if Gene_dict[circuit_plasmid][item]['feature'] in ['gene','resistance_maker']:
						start = Gene_dict[circuit_plasmid][item]['start']
						end = Gene_dict[circuit_plasmid][item]['end']
						if Gene_dict[circuit_plasmid][item]['strand'] == '+':
							circuit_gene_dict_P.append({'NAME': item, 'START': start, 'END': end})
						elif Gene_dict[circuit_plasmid][item]['strand'] == '-':
							circuit_gene_dict_M.append({'NAME': item, 'START': start, 'END': end})

				reporter_gene_dict_P = []
				reporter_gene_dict_M = []
				for item in Gene_dict[reporter_plasmid]:
					if Gene_dict[reporter_plasmid][item]['feature'] in ['gene','resistance_maker']:
						start = Gene_dict[reporter_plasmid][item]['start']
						end = Gene_dict[reporter_plasmid][item]['end']
						if Gene_dict[reporter_plasmid][item]['strand'] == '+':
							reporter_gene_dict_P.append({'NAME': item, 'START': start, 'END': end})
						elif Gene_dict[reporter_plasmid][item]['strand'] == '-':
							reporter_gene_dict_M.append({'NAME': item, 'START': start, 'END': end})
				
				genome_gene_dict_P = []
				genome_gene_dict_M = []
				for item in Gene_dict[genome]:
					if Gene_dict[genome][item]['feature'] in ['gene','resistance_maker']:
						start = Gene_dict[genome][item]['start']
						end = Gene_dict[genome][item]['end']
						if Gene_dict[genome][item]['strand'] == '+':
							genome_gene_dict_P.append({'NAME': item, 'START': start, 'END': end})
						elif Gene_dict[genome][item]['strand'] == '-':
							genome_gene_dict_M.append({'NAME': item, 'START': start, 'END': end})
				
				# ----------------------------------- #
				
				circuit_reads_dict_P = Ribo_seq_Profile[sample_name][circuit_plasmid]['fwd']
				circuit_reads_dict_M = Ribo_seq_Profile[sample_name][circuit_plasmid]['rev']
				
				reporter_reads_dict_P = Ribo_seq_Profile[sample_name][reporter_plasmid]['fwd']
				reporter_reads_dict_M = Ribo_seq_Profile[sample_name][reporter_plasmid]['rev']
				
				genome_reads_dict_P = Ribo_seq_Profile[sample_name][genome]['fwd']
				genome_reads_dict_M = Ribo_seq_Profile[sample_name][genome]['rev']
				
				circuit_dropoff_corrected_reads_dict_P, circuit_dropoff_corrected_reads_dict_M = self.exp_dropoff_correction(genome_gene_dict_P, genome_gene_dict_M, genome_reads_dict_P, genome_reads_dict_M, circuit_gene_dict_P, circuit_gene_dict_M, circuit_reads_dict_P, circuit_reads_dict_M, 1.0*(self.total_active_ribosomes/1e9), 150, 5, [.446, 6e-3, 1], self.circuit_output_RD_exponential_fit_file[sample_name])
				circuit_pause_corrected_reads_dict_P, circuit_pause_corrected_reads_dict_M = self.pause_site_correction(circuit_gene_dict_P, circuit_gene_dict_M, circuit_dropoff_corrected_reads_dict_P, circuit_dropoff_corrected_reads_dict_M)
			
				reporter_dropoff_corrected_reads_dict_P, reporter_dropoff_corrected_reads_dict_M = self.exp_dropoff_correction(genome_gene_dict_P, genome_gene_dict_M, genome_reads_dict_P, genome_reads_dict_M, reporter_gene_dict_P, reporter_gene_dict_M, reporter_reads_dict_P, reporter_reads_dict_M, 1.0*(self.total_active_ribosomes/1e9), 150, 5, [.446, 6e-3, 1], self.reporter_output_RD_exponential_fit_file[sample_name])
				reporter_pause_corrected_reads_dict_P, reporter_pause_corrected_reads_dict_M = self.pause_site_correction(reporter_gene_dict_P, reporter_gene_dict_M, reporter_dropoff_corrected_reads_dict_P, reporter_dropoff_corrected_reads_dict_M)
				
				genome_dropoff_corrected_reads_dict_P, genome_dropoff_corrected_reads_dict_M = self.exp_dropoff_correction(genome_gene_dict_P, genome_gene_dict_M, genome_reads_dict_P, genome_reads_dict_M, genome_gene_dict_P, genome_gene_dict_M, genome_reads_dict_P, genome_reads_dict_M, 1.0*(self.total_active_ribosomes/1e9), 150, 5, [.446, 6e-3, 1], self.genome_output_RD_exponential_fit_file[sample_name])
				genome_pause_corrected_reads_dict_P, genome_pause_corrected_reads_dict_M = self.pause_site_correction(genome_gene_dict_P, genome_gene_dict_M, genome_dropoff_corrected_reads_dict_P, genome_dropoff_corrected_reads_dict_M)
				
				# ----------------------------------- #
				
				for gene in circuit_gene_dict_P: 
					start = gene['START'] + aa_ends_excluded*3 - 1 #excludes first and last 5 codons to remove effects of translation initiaon and termination
					stop  = gene['END'] - aa_ends_excluded*3
					RD_dict[sample_name][circuit_plasmid][gene['NAME']] = np.mean(circuit_pause_corrected_reads_dict_P[start:stop+1])
					
				for gene in circuit_gene_dict_M: 
					start = gene['START'] + aa_ends_excluded*3 - 1 #excludes first and last 5 codons to remove effects of translation initiaon and termination
					stop  = gene['END'] - aa_ends_excluded*3
					RD_dict[sample_name][circuit_plasmid][gene['NAME']] = np.mean(circuit_pause_corrected_reads_dict_M[start:stop+1])
					
				for gene in reporter_gene_dict_P: 
					start = gene['START'] + aa_ends_excluded*3 - 1 #excludes first and last 5 codons to remove effects of translation initiaon and termination
					stop  = gene['END'] - aa_ends_excluded*3
					RD_dict[sample_name][reporter_plasmid][gene['NAME']] = np.mean(reporter_pause_corrected_reads_dict_P[start:stop+1])
					
				for gene in reporter_gene_dict_M: 
					start = gene['START'] + aa_ends_excluded*3 - 1 #excludes first and last 5 codons to remove effects of translation initiaon and termination
					stop  = gene['END'] - aa_ends_excluded*3
					RD_dict[sample_name][reporter_plasmid][gene['NAME']] = np.mean(reporter_pause_corrected_reads_dict_M[start:stop+1])
					
				for gene in genome_gene_dict_P: 
					start = gene['START'] + aa_ends_excluded*3 - 1 #excludes first and last 5 codons to remove effects of translation initiaon and termination
					stop  = gene['END'] - aa_ends_excluded*3
					RD_dict[sample_name][genome][gene['NAME']] = np.mean(genome_pause_corrected_reads_dict_P[start:stop+1])
					
				for gene in genome_gene_dict_M: 
					start = gene['START'] + aa_ends_excluded*3 - 1 #excludes first and last 5 codons to remove effects of translation initiaon and termination
					stop  = gene['END'] - aa_ends_excluded*3
					RD_dict[sample_name][genome][gene['NAME']] = np.mean(genome_pause_corrected_reads_dict_M[start:stop+1])
				
		# ----------------------------------- #
		
		workbook = xlwt.Workbook()
		
		worksheet = workbook.add_sheet('RD')
		
		worksheet.write(0,0,"Reference")
		worksheet.write(0,1,"Gene")
		
		[circuit_plasmid, reporter_plasmid, genome] = self.reference_name['RiboSeq_State_1']
		
		row2 = 1
		for chromos in self.reference_name['RiboSeq_State_1']:
			for gene2 in sorted(Gene_dict[chromos]):
				if Gene_dict[chromos][gene2]['feature'] in ['gene','resistance_maker']:
					if chromos == circuit_plasmid:
						worksheet.write(row2,0,'Circuit_Plasmid')
					elif chromos == reporter_plasmid:
						worksheet.write(row2,0,'Reporter_Plasmid')
					elif chromos == genome:
						worksheet.write(row2,0,'NC_010473.1')
					worksheet.write(row2,1,gene2.split()[0])
					row2 += 1
		
		kk = 2
		for sample_name in sorted(RD_dict):
			worksheet.write(0,kk,sample_name)
			row2 = 1
			for chromos in self.reference_name['RiboSeq_State_1']:
				for gene2 in sorted(Gene_dict[chromos]):
					if Gene_dict[chromos][gene2]['feature'] in ['gene','resistance_maker']:
						if sample_name == 'RiboSeq_Control':
							[circuit_plasmid2, reporter_plasmid2, genome2] = self.reference_name[sample_name]
							
							if chromos == circuit_plasmid:
								chromos2 = circuit_plasmid2
							elif chromos == reporter_plasmid:
								chromos2 = reporter_plasmid2
							else:
								chromos2 = chromos
							
							if gene2 in RD_dict[sample_name][chromos2]:
								worksheet.write(row2,kk,RD_dict[sample_name][chromos2][gene2])
							else:
								worksheet.write(row2,kk,0)
							row2 += 1
						else:
							worksheet.write(row2,kk,RD_dict[sample_name][chromos][gene2])
							row2 += 1
			
			kk += 1
			
		workbook.save(self.output_RD_table_filename+'/RD.xls')
		
		return RD_dict
	
	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def gene_length_dropoff_function(self, gene_dict_P, gene_dict_M, reads_dict_P, reads_dict_M, min_read_count, min_length, aa_ends_excluded, p0, output_expo):
		
		"""Function: Tries to check for dropoff of read densities along gene by comparing read densities in 
		50 codon windows against first 50 codons. This function tries to fit the dropoff of read densities with 
		an exponential curve, plots the read densities against the fitted curve and returns the parameters 
		Input1(genes): gene dictionary of genes with their start and end position
		Input2(reads_dict): read densities along position on the genome
		Input3(min_read_count): Float. The minimum number of reads to qualify a gene for inclusion. Default = 1*(self.total_active_ribosomes/1e9)
		Input4(min_length): Integer. Minimum length of gene to qualify for inclusion. Default = 150
		Input5(aa_ends_excluded): Number of amino acids to exclude from each end. Default = 5
		Input6(p0): List of 3 floats. Input value for curve fitting optimization. Defulat value = [0.446, 6e-3, 1], 
		This is specific for E coli. 
		Output1,2,and3 (A, K, C): float values for fitted parameters to the dropoff curve. 
		"""
		
		window_counter = {1:[]}   # all genes we look at will have first window
		window_averages = []
		
		for gene in gene_dict_P:
			
			start = gene['START'] + aa_ends_excluded*3 - 1 #excludes first and last 5 codons 
			stop  = gene['END'] - aa_ends_excluded*3
			length = stop - start + 1
			gene_pos = list(range(start, stop+1))
			
			#readdensity = sum(reads_dict_P[start:stop+1])
			readdensity = np.mean(reads_dict_P[start:stop+1])
			
			if length > min_length and readdensity > min_read_count:
				
				window = 1   # starting from second window 
				pos_in_win = list(range(start, start + min_length + 1))
				
				firstwindow = sum(reads_dict_P[start:start + min_length + 1])
				
				if firstwindow == 0.:
					continue

				window_counter[window].append(1.)
				window += 1
				while start + min_length * (window) < stop: 
					
					pos_in_ea_win = list(range(start + min_length * (window-1), start + min_length * (window) + 1))
					
					reads_in_ea_win = sum(reads_dict_P[start + min_length * (window-1):start + min_length * (window) + 1])
					
					if window in window_counter:
						window_counter[window].append(reads_in_ea_win/firstwindow)
					else: 
						window_counter[window] = [reads_in_ea_win/firstwindow]
					
					window += 1
					
						
		for gene in gene_dict_M:
			
			start = gene['START'] + aa_ends_excluded*3 - 1 #excludes first and last 5 codons 
			stop  = gene['END'] - aa_ends_excluded*3
			length = stop - start + 1
			gene_pos = list(range(start, stop+1))
			
			#readdensity = sum(reads_dict_M[start:stop+1])
			readdensity = np.mean(reads_dict_M[start:stop+1])
			
			if length > min_length and readdensity > min_read_count:
				
				window = 1 
				pos_in_win = list(range(stop - min_length, stop+1))
				
				firstwindow = sum(reads_dict_M[stop - min_length:stop+1])
				
				if firstwindow == 0.:
					continue
				
				window_counter[window].append(1.)
				window += 1
				while stop - min_length * (window) > stop: 
					
					pos_in_ea_win = list(range(stop - min_length * (window + 1), stop - min_length * (window) + 1))
					
					reads_in_ea_win = sum(reads_dict_M[stop - min_length * (window + 1):stop - min_length * (window) + 1])
					
					if window in window_counter: 
						window_counter[window].append(reads_in_ea_win/firstwindow)
					else: 
						window_counter[window] = [reads_in_ea_win/firstwindow]
						
					window += 1


		for i in range(1, len(window_counter) + 1):
			window_averages.append(np.median(window_counter[i]))

		plotx = range(aa_ends_excluded*3, 875, min_length)
		plot_x = np.array(plotx)
		popt, pcov = optimize.curve_fit(self.model_func, plot_x, window_averages[0:len(plotx)], p0=p0)
		A, K, C = popt[0], popt[1], popt[2]
		
		residuals = window_averages[0:len(plotx)] - self.model_func(plot_x, A, K, C)
		sqr_res = sum(residuals**2)
		print('The calculated function parameters are: \nA=%f\nK=%f\nC=%f\nSum of residuals squared=%f'%(A,K,C,sqr_res))
		
		fig = plt.figure()
		ax = fig.add_subplot(111)
		plt.scatter(plotx, window_averages[0:len(plotx)], color='forestgreen')
		plt.plot(plot_x, self.model_func(plot_x, A, K, C), color='firebrick')
		plt.ylim([0, 1.1])
		plt.xlabel('Distance from start codon (codons)')
		plt.ylabel('Relative ribosome occupancy')
		plt.show()
		fig.savefig(output_expo)
		
		return A, K, C

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def exp_dropoff_correction(self, endo_gene_dict_P, endo_gene_dict_M, endo_reads_dict_P, endo_reads_dict_M, gene_dict_P, gene_dict_M, reads_dict_P, reads_dict_M, min_read_count, min_length, aa_ends_excluded, p0, outputfile):
		
		A, K, C = self.gene_length_dropoff_function(endo_gene_dict_P, endo_gene_dict_M, endo_reads_dict_P, endo_reads_dict_M, min_read_count, min_length, aa_ends_excluded, p0, outputfile)

		i = range(1, 10000)
		multiplier = []
		for nt in i:
			j = self.model_func(float(nt), A, K, C)
			multiplier.append(j)

		dropoff_corrected_reads_dict_P = reads_dict_P.copy() 
		dropoff_corrected_reads_dict_M = reads_dict_M.copy() 

		# plus strand 
		for gene in gene_dict_P: 
			
			start = gene['START'] + aa_ends_excluded*3 - 1 #exclude first and last 5 codons
			stop  = gene['END'] - aa_ends_excluded*3
			length = stop - start + 1

			for i in range(length):
				dropoff_corrected_reads_dict_P[start + i] = dropoff_corrected_reads_dict_P[start + i] / multiplier[i]

		# minus strand 
		for gene in gene_dict_M: 
			
			start = gene['START'] + aa_ends_excluded*3 - 1
			stop  = gene['END'] - aa_ends_excluded*3
			length = stop - start + 1

			for i in range(length):
				dropoff_corrected_reads_dict_M[stop - i] = dropoff_corrected_reads_dict_M[stop - i] / multiplier[i]

		return dropoff_corrected_reads_dict_P, dropoff_corrected_reads_dict_M 

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def pause_site_correction(self, gene_dict_P, gene_dict_M, reads_dict_P, reads_dict_M, method = 'std_dev', set_to_zero = False, aa_ends_excluded = 5):
		
		pause_corrected_reads_dict_P = reads_dict_P.copy() 
		pause_corrected_reads_dict_M = reads_dict_M.copy()
		
		# correct for plus strand 
		for i in range(2):
			
			gene_dict  = [gene_dict_P, gene_dict_M][i]
			reads_dict = [reads_dict_P, reads_dict_M][i]
			pause_corrected_reads_dict = [pause_corrected_reads_dict_P, pause_corrected_reads_dict_M][i]
			
			for gene in gene_dict: 
				
				start = gene['START'] + aa_ends_excluded*3 - 1
				stop  = gene['END'] - aa_ends_excluded*3 
				length = stop - start + 1
				gene_pos = list(range(start, stop+1))
				
				readdensity = reads_dict[start:stop+1])

				avg = sum(readdensity)/length # find the average occupancy across the gene 
				std_dev = np.std(readdensity)  # find the std
				
				if method == 'std_dev':
					threshold = avg + 5*std_dev
					
				for i in range(length):
					
					try:
						if reads_dict[start+i] < threshold: #pause sites are considered to be nucleotides where occupancy is > threshold 
							pause_corrected_reads_dict[start+i] = reads_dict[start+i]
							
						elif set_to_zero: # occupany on identified pause sites are reduced to threshold value
							pause_corrected_reads_dict[start+i] = 0
							
						else: 
							pause_corrected_reads_dict[start+i] = threshold
							
					except KeyError: 
						print(gene)
						
		return pause_corrected_reads_dict_P, pause_corrected_reads_dict_M

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def calculate_MW(self, FASTA_filename, GFF_filename):
		
		Genome_Sequence = {}
		with open(FASTA_filename, 'rU') as ins:
			for lines in ins:
				field=lines.strip().split()
				#print field
				if field != []:
					if field[0].count('>') == 1:
						ID = field[0].split('>')[1]
						Genome_Sequence[ID] = []
				else:
					Genome_Sequence[ID].append(field[0])
		
		# ----------------------------------- #
		
		for item in Genome_Sequence:
			Genome_Sequence[item] = ''.join(Genome_Sequence[item]).upper()
		
		# ----------------------------------- #
		
		Genome_CDS = {item:{} for item in Genome_Sequence}
		with open(GFF_filename, 'rU') as ins:
			for lines in ins:
				field=lines.strip().split()
				if field != []:
					if field[2] == 'gene' or field[2] == 'resistance_maker':
						genome, start, end, strand, CDS_name = field[0], field[3],field[4],field[6],field[8].split('=')[1]
						Genome_CDS[genome][CDS_name] = [start, end, strand]
		#print Genome_CDS
		
		# ----------------------------------- #
		
		CDS_seq = {item:{} for item in Genome_Sequence}
		AA_seq = {item:{} for item in Genome_Sequence}
		for genome in Genome_CDS:
			for xxx in Genome_CDS[genome]:
				#print xxx,Genome_CDS[genome][xxx]
				if Genome_CDS[genome][xxx][2] == '+':
					CDS_seq[genome][xxx] = Genome_Sequence[genome][int(Genome_CDS[genome][xxx][0])-1 : int(Genome_CDS[genome][xxx][1])].upper()
					AA_seq[genome][xxx] = self.amino_acid_seq (CDS_seq[genome][xxx], codon_table)
				elif Genome_CDS[genome][xxx][2] == '-':
					temp_seq = Genome_Sequence[genome][int(Genome_CDS[genome][xxx][0])-1 : int(Genome_CDS[genome][xxx][1])].upper()
					CDS_seq[genome][xxx] = self.reverse_complement_naive (temp_seq)
					AA_seq[genome][xxx] = self.amino_acid_seq (CDS_seq[genome][xxx], codon_table)
		#print CDS_seq
		#print AA_seq
		
		# ----------------------------------- #
		
		amino_acid_usage = {item:{} for item in Genome_Sequence}
		for genome in AA_seq:
			for gene in AA_seq[genome]:
				amino_acid_usage[genome][gene] = {}
				for key in codon_table:
					if key != 'STOP':
						amino_acid_usage[genome][gene][key] = AA_seq[genome][gene].count(key)
		#print amino_acid_usage
		
		# ----------------------------------- #
		
		Protein_MW = {item:{} for item in Genome_Sequence}
		for genome in amino_acid_usage:
			for gene in amino_acid_usage[genome]:
				Protein_MW[genome][gene] = 0
				for key in amino_acid_usage[genome][gene]:
					Protein_MW[genome][gene] += amino_acid_usage[genome][gene][key]*amino_acid_MW[key]
				
				Protein_MW[genome][gene] -= (sum(amino_acid_usage[genome][gene].values())-1)*18.01528		# to account for loss of water for evey peptide bond formation
		#print Protein_MW
		
		return Protein_MW

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def generate_MW_table(self, Protein_MW):
		
		workbook = xlwt.Workbook()
		
		worksheet = workbook.add_sheet('Molecular_Weight')

		worksheet.write(0,0,"Reference")
		worksheet.write(0,1,"Gene")
		worksheet.write(0,2,"Molecular_weight_(g/mol)")
		
		row2 = 1
		
		for chromos in Protein_MW:
			for gene2 in sorted(Protein_MW[chromos]):
				if chromos == '0x41v70':
					worksheet.write(row2,0,'Circuit_Plasmid')
				elif chromos == 'Dv_pBM3R1-YFP':
					worksheet.write(row2,0,'Reporter_Plasmid')
				elif chromos == 'NC_010473.1':
					worksheet.write(row2,0,'NC_010473.1')
				worksheet.write(row2,1,gene2)
				worksheet.write(row2,2,Protein_MW[chromos][gene2])
				row2 += 1
				
		workbook.save(self.output_Protein_MW_table_filename+'/Molecular_Weight.xls')
				
	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def calculate_proteome_fraction(self, RD, Protein_MW):

		Total_proteome = {}
		Proteome = {}
		for sample_name in RD:
			temp = 0.0
			Proteome[sample_name] = {}
			for genome in RD[sample_name]:
				Proteome[sample_name][genome] = {}
				for gene in RD[sample_name][genome]:
					Proteome[sample_name][genome][gene] = RD[sample_name][genome][gene]*Protein_MW[genome][gene]
					temp += RD[sample_name][genome][gene]*Protein_MW[genome][gene]
			Total_proteome[sample_name] = temp
		
		Proteome_Fraction = {}
		for sample_name in Proteome:
			Proteome_Fraction[sample_name] = {}
			for genome in Proteome[sample_name]:
				Proteome_Fraction[sample_name][genome] = {}
				for gene in Proteome[sample_name][genome]:
					Proteome_Fraction[sample_name][genome][gene] = Proteome[sample_name][genome][gene]/Total_proteome[sample_name]

		return Proteome_Fraction

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def generate_proteome_fraction_table(self, Proteome_Fraction):

		workbook = xlwt.Workbook()

		worksheet = workbook.add_sheet('Proteome_Fraction')

		worksheet.write(0,0,"Reference")
		worksheet.write(0,1,"Gene")

		[circuit_plasmid, reporter_plasmid, genome] = self.reference_name['RiboSeq_State_1']

		row2 = 1
		for chromos in self.reference_name['RiboSeq_State_1']:
			for gene2 in sorted(Proteome_Fraction['RiboSeq_State_1'][chromos]):
				if chromos == circuit_plasmid:
					worksheet.write(row2,0,'Circuit_Plasmid')
				elif chromos == reporter_plasmid:
					worksheet.write(row2,0,'Reporter_Plasmid')
				elif chromos == genome:
					worksheet.write(row2,0,'NC_010473.1')
				worksheet.write(row2,1,gene2)
				row2 += 1
		
		kk = 2
		for sample_name in sorted(Proteome_Fraction):
			worksheet.write(0,kk,sample_name)
			row2 = 1
			for chromos in self.reference_name['RiboSeq_State_1']:
				for gene2 in sorted(Proteome_Fraction['RiboSeq_State_1'][chromos]):
					worksheet.write(row2,kk,Proteome_Fraction[sample_name][chromos][gene2])
					row2 += 1
			
			kk += 1
			
		workbook.save(self.output_Proteome_Fraction_table_filename+'/Proteome_Fraction.xls')

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def generate_gate_response_function(self, RNA_seq_Profile, Ribo_seq_Profile, alfa, alfa2, promoter_TSS, promoter_TSS2, X2, X1, Gates, RBSs, terminators, gate_list):
		
		gamma = self.RNAP_p15A_total		# to convert to RNAP / s
		gamma2 = self.RNAP_p15A_total		# to convert to RNAP / s
		
		Gate_input_Cello = {}
		Gate_output_Cello = {}
		Gate_input_RNAP_Flux = {}
		Gate_output_RNAP_Flux = {}
		Gate_ribosome_usage = {}

		for gtg in range(len(gate_list)):
			
			#-------------------------------------------------------------------------------------------------------------#
			
			predit_flux = np.logspace(-6,np.log10(75),1000) * gamma
			
			if gate_list[gtg] == 'BM3R1':
				Cello = 0.4 * (self.par['p'+gate_list[gtg]]['min']*gamma2+(self.par['p'+gate_list[gtg]]['max']*gamma2-self.par['p'+gate_list[gtg]]['min']*gamma2)*(((self.par['p'+gate_list[gtg]]['K']*gamma2)**self.par['p'+gate_list[gtg]]['n'])/(((self.par['p'+gate_list[gtg]]['K']*gamma2)**self.par['p'+gate_list[gtg]]['n'])+(predit_flux**self.par['p'+gate_list[gtg]]['n']))) )
			else:
				Cello = self.par['p'+gate_list[gtg]]['min']*gamma+(self.par['p'+gate_list[gtg]]['max']*gamma-self.par['p'+gate_list[gtg]]['min']*gamma)*(((self.par['p'+gate_list[gtg]]['K']*gamma)**self.par['p'+gate_list[gtg]]['n'])/(((self.par['p'+gate_list[gtg]]['K']*gamma)**self.par['p'+gate_list[gtg]]['n'])+(predit_flux**self.par['p'+gate_list[gtg]]['n'])))
				
			#-------------------------------------------------------------------------------------------------------------#
			
			input_list_avg = []
			output_list_avg = []
			ribosome_list_avg = []

			for uu in range(self.inducer):
				print "-------------state-------------", uu
				
				[circuit_plasmid, reporter_plasmid, genome] = self.reference_name['RNASeq_State_'+str(uu+1)]
				
				flux_dic_avg = {}
				for item_cir in promoter_TSS:
					
					profile = RNA_seq_Profile['RNASeq_State_'+str(uu+1)][circuit_plasmid]['fwd'][promoter_TSS[item_cir]-X2:promoter_TSS[item_cir]-X1]
					flux_dic_before_avg = np.mean([alfa * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])
					
					profile = RNA_seq_Profile['RNASeq_State_'+str(uu+1)][circuit_plasmid]['fwd'][promoter_TSS[item_cir]+X1:promoter_TSS[item_cir]+X2]
					flux_dic_after_avg = np.mean([alfa * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])
									
					flux_dic_avg[item_cir] = max(1e-6, flux_dic_after_avg - flux_dic_before_avg)
					
				#-------------------------------------------------------------------------------------------------------------#
					
				flux_dic_avg_2 = {}
				for item_rep in promoter_TSS2:

					profile = RNA_seq_Profile['RNASeq_State_'+str(uu+1)][reporter_plasmid]['fwd'][promoter_TSS2[item_rep]-X2:promoter_TSS2[item_rep]-X1]
					flux_dic_before_avg_2 = np.mean([alfa2 * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])

					profile = RNA_seq_Profile['RNASeq_State_'+str(uu+1)][reporter_plasmid]['fwd'][promoter_TSS2[item_rep]+X1:promoter_TSS2[item_rep]+X2]
					flux_dic_after_avg_2 = np.mean([alfa2 * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])

					flux_dic_avg_2[item_rep] = max(1e-6, flux_dic_after_avg_2 - flux_dic_before_avg_2)
					
				#-------------------------------------------------------------------------------------------------------------#

				input_activity_avg = sum([flux_dic_avg[fgh] for fgh in Gates[gate_list[gtg]]])

				if gate_list[gtg] == 'BM3R1':
					output_activity_avg = flux_dic_avg_2['pBM3R1']
				else:
					output_activity_avg = flux_dic_avg['p'+gate_list[gtg]]

				input_list_avg.append(input_activity_avg)
				output_list_avg.append(output_activity_avg)
				
				#-------------------------------------------------------------------------------------------------------------#
				
				avg_Ribo = np.sum(np.array(Ribo_seq_Profile['RiboSeq_State_'+str(uu+1)][circuit_plasmid][RBSs[gtg][0]:terminators[gtg][0]]))
				
				ribosome_list_avg.append(avg_Ribo)
				
				#-------------------------------------------------------------------------------------------------------------#
				
			Gate_input_Cello[gate_list[gtg]] = predit_flux
			Gate_output_Cello[gate_list[gtg]] = Cello
			Gate_input_RNAP_Flux[gate_list[gtg]] = input_list_avg
			Gate_output_RNAP_Flux[gate_list[gtg]] = output_list_avg
			Gate_ribosome_usage[gate_list[gtg]] = ribosome_list_avg
			
		return Gate_input_Cello, Gate_output_Cello, Gate_input_RNAP_Flux, Gate_output_RNAP_Flux, Gate_ribosome_usage
			
	# ---------------------------------------------------------------------------------------------------------------------------------------- #
		
	def calculate_circuit_total_RNAP_flux(self, RNA_seq_Profile, alfa, alfa2, border_circuit_control, border_reporter_control, border_circuit, border_reporter):
		
		Total_RNAP_Flux = {}
		
		for sample_name in RNA_seq_Profile:
			
			[circuit_plasmid, reporter_plasmid, genome] = self.reference_name[sample_name]
			
			if sample_name == 'RNASeq_Control':
				
				profile = RNA_seq_Profile[sample_name][circuit_plasmid]['fwd'][border_circuit_control[0]:border_circuit_control[1]]
				RNAP_FLux_circuit_fwd = sum([alfa * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])
				
				profile = RNA_seq_Profile[sample_name][circuit_plasmid]['rev'][border_circuit_control[0]:border_circuit_control[1]]
				RNAP_FLux_circuit_rev = sum([alfa * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])
				
				profile = RNA_seq_Profile[sample_name][reporter_plasmid]['fwd'][border_reporter_control[0]:border_reporter_control[1]]
				RNAP_FLux_reporter_fwd = sum([alfa2 * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])
				
				profile = RNA_seq_Profile[sample_name][reporter_plasmid]['rev'][border_reporter_control[0]:border_reporter_control[1]]
				RNAP_FLux_reporter_rev = sum([alfa2 * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])
				
			else:
					
				profile = RNA_seq_Profile[sample_name][circuit_plasmid]['fwd'][border_circuit[0]:border_circuit[1]]
				RNAP_FLux_circuit_fwd = sum([alfa * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])
				
				profile = RNA_seq_Profile[sample_name][circuit_plasmid]['rev'][border_circuit[0]:border_circuit[1]]
				RNAP_FLux_circuit_rev = sum([alfa * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])
				
				profile = RNA_seq_Profile[sample_name][reporter_plasmid]['fwd'][border_reporter[0]:border_reporter[1]]
				RNAP_FLux_reporter_fwd = sum([alfa2 * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])
				
				profile = RNA_seq_Profile[sample_name][reporter_plasmid]['rev'][border_reporter[0]:border_reporter[1]]
				RNAP_FLux_reporter_rev = sum([alfa2 * self.RNAP_p15A_total * (10.0**(self.coef1 * math.log(self.degradation * x, 10) + self.coef2)) for x in profile])
				
			Total_RNAP_Flux[sample_name] = RNAP_FLux_circuit_fwd + RNAP_FLux_circuit_rev + RNAP_FLux_reporter_fwd + RNAP_FLux_reporter_rev

		return Total_RNAP_Flux

	# ---------------------------------------------------------------------------------------------------------------------------------------- #

	def calculate_total_ribosome_usage(self, Ribo_seq_Profile, border_circuit_control, border_reporter_control, border_circuit, border_reporter):

		Total_Ribo_usage = {}

		for sample_name in Ribo_seq_Profile:
			
			[circuit_plasmid, reporter_plasmid, genome] = self.reference_name[sample_name]
			
			if sample_name == 'RiboSeq_Control':

				Ribo_usage_circuit_fwd = np.sum(np.array(Ribo_seq_Profile[sample_name][circuit_plasmid]['fwd'][border_circuit_control[0]:border_circuit_control[1]]))

				Ribo_usage_circuit_rev = np.sum(np.array(Ribo_seq_Profile[sample_name][circuit_plasmid]['rev'][border_circuit_control[0]:border_circuit_control[1]]))

				Ribo_usage_reporter_fwd = np.sum(np.array(Ribo_seq_Profile[sample_name][reporter_plasmid]['fwd'][border_reporter_control[0]:border_reporter_control[1]]))

				Ribo_usage_reporter_rev = np.sum(np.array(Ribo_seq_Profile[sample_name][reporter_plasmid]['rev'][border_reporter_control[0]:border_reporter_control[1]]))

			else:
				
				Ribo_usage_circuit_fwd = np.sum(np.array(Ribo_seq_Profile[sample_name][circuit_plasmid]['fwd'][border_circuit[0]:border_circuit[1]]))

				Ribo_usage_circuit_rev = np.sum(np.array(Ribo_seq_Profile[sample_name][circuit_plasmid]['rev'][border_circuit[0]:border_circuit[1]]))

				Ribo_usage_reporter_fwd = np.sum(np.array(Ribo_seq_Profile[sample_name][reporter_plasmid]['fwd'][border_reporter[0]:border_reporter[1]]))

				Ribo_usage_reporter_rev = np.sum(np.array(Ribo_seq_Profile[sample_name][reporter_plasmid]['rev'][border_reporter[0]:border_reporter[1]]))

			Total_Ribo_usage[sample_name] = Ribo_usage_circuit_fwd + Ribo_usage_circuit_rev + Ribo_usage_reporter_fwd + Ribo_usage_reporter_rev

		return Total_Ribo_usage

# ---------------------------------------------------------------------------------------------------------------------------------------- #

if __name__ == "__main__":
	
	# ---------------------------------- #
	
	par = {}

	""" from Table S6 in Cello Paper """
	par['pTac']  = {'off':0.0034,'on':2.8}
	par['pTet1'] = {'off':0.0013,'on':4.4}
	par['pTet2'] = {'off':0.0013,'on':4.4}
	par['pBAD1'] = {'off':0.0082,'on':2.5}
	par['pBAD2'] = {'off':0.0082,'on':2.5}

	""" from Table S4 in Cello Paper """
	par['pAmtR']    = {'min':0.06, 'max':3.8,'K':0.07,'n':1.6}
	par['pBM3R1']   = {'min':0.005,'max':0.5,'K':0.15,'n':2.9}	# with RBS B2
	par['pSrpR']    = {'min':0.007,'max':2.1,'K':0.10,'n':2.8}	# with RBS S4
	par['pPhlF']    = {'min':0.02, 'max':4.1,'K':0.13,'n':3.9}	# with RBS P2
	par['pBetI']    = {'min':0.07, 'max':3.8,'K':0.41,'n':2.4}
	par['pHlyIIR']  = {'min':0.07, 'max':2.5,'K':0.19,'n':2.6}
	par['pAmeR']    = {'min':0.20, 'max':3.8,'K':0.09,'n':1.4}

	# ---------------------------------- #

	inducer = {}
	inducer[1] = {'pTac':'off','pTet1':'off','pBAD1':'off','pTet2':'off','pBAD2':'off'}
	inducer[2] = {'pTac':'on' ,'pTet1':'off','pBAD1':'off','pTet2':'off','pBAD2':'off'}
	inducer[3] = {'pTac':'off','pTet1':'on' ,'pBAD1':'off','pTet2':'on' ,'pBAD2':'off'}
	inducer[4] = {'pTac':'on' ,'pTet1':'on' ,'pBAD1':'off','pTet2':'on' ,'pBAD2':'off'}
	inducer[5] = {'pTac':'off','pTet1':'off','pBAD1':'on' ,'pTet2':'off','pBAD2':'on' }
	inducer[6] = {'pTac':'on' ,'pTet1':'off','pBAD1':'on' ,'pTet2':'off','pBAD2':'on' }
	inducer[7] = {'pTac':'off','pTet1':'on' ,'pBAD1':'on' ,'pTet2':'on' ,'pBAD2':'on' }
	inducer[8] = {'pTac':'on' ,'pTet1':'on' ,'pBAD1':'on' ,'pTet2':'on' ,'pBAD2':'on' }
	
	# ---------------------------------- #
	
	degradation = 0.0067
	
	RNAP_p15A_total = 0.171	#RNAP/s   (on p15A)
	RNAP_p15A = 0.019		#RNAP/s-DNA
	
	coef1 = 1.6369
	coef2 = -4.2971
	
	omega = 15			# codon / s	ribosome translocation rate
	
	total_active_ribosomes = 20000
	
	# ---------------------------------- #
		
	codon_table = {'phe': ['TTT','TTC'],
			'leu': ['TTA','TTG','CTT','CTC','CTA','CTG'],
			'ile': ['ATT','ATC','ATA'],
			'met': ['ATG'],
			'val': ['GTT','GTC','GTA','GTG'],
			'ser': ['TCT','TCC','TCA','TCG','AGT','AGC'],
			'pro': ['CCT','CCC','CCA','CCG'],
			'thr': ['ACT','ACC','ACA','ACG'],
			'ala': ['GCT','GCC','GCA','GCG'],
			'tyr': ['TAT','TAC'],
			'STOP': ['TAA','TAG','TGA'],
			'his': ['CAT','CAC'],
			'gln': ['CAA','CAG'],
			'asn': ['AAT','AAC'],
			'lys': ['AAA','AAG'],
			'asp': ['GAT','GAC'],
			'glt': ['GAA','GAG'],
			'cys': ['TGT','TGC'],
			'trp': ['TGG'],
			'arg': ['CGT','CGC','CGA','CGG','AGA','AGG'],
			'gly': ['GGT','GGC','GGA','GGG']}
		
	# ---------------------------------- #
		
	amino_acid_MW = {'ala':89.094,
			'arg':174.204,
			'asn':132.119,
			'asp':133.103,
			'cys':121.15,
			'gln':146.146,
			'glt':147.130,
			'gly':75.067,
			'his':155.157,
			'ile':131.175,
			'leu':131.175,
			'lys':146.190,
			'met':149.21,
			'phe':165.192,
			'pro':115.132,
			'ser':105.093,
			'thr':119.120,
			'trp':204.229,
			'tyr':181.191,
			'val':117.148}
		
	# ---------------------------------- #
	
	frag_dist_bin = (0,50)
	
	# ---------------------------------- #
	
	all_sample_names =  (
	'RNASeq_Control','RNASeq_State_1','RNASeq_State_2','RNASeq_State_3','RNASeq_State_4','RNASeq_State_5','RNASeq_State_6','RNASeq_State_7','RNASeq_State_8',
	'RiboSeq_Control','RiboSeq_State_1','RiboSeq_State_2','RiboSeq_State_3','RiboSeq_State_4','RiboSeq_State_5','RiboSeq_State_6','RiboSeq_State_7','RiboSeq_State_8'
	)
	
	# ---------------------------------- #
	
	reference_name = { 
			'RNASeq_Control':('pAN871', 'Dv4_pBAD-YFP', 'NC_010473.1'),
			'RNASeq_State_1':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RNASeq_State_2':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RNASeq_State_3':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RNASeq_State_4':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RNASeq_State_5':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RNASeq_State_6':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RNASeq_State_7':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RNASeq_State_8':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RiboSeq_Control':('pAN871', 'Dv4_pBAD-YFP', 'NC_010473.1'),
			'RiboSeq_State_1':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RiboSeq_State_2':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RiboSeq_State_3':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RiboSeq_State_4':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RiboSeq_State_5':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RiboSeq_State_6':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RiboSeq_State_7':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1'),
			'RiboSeq_State_8':('0x41v70', 'Dv_pBM3R1-YFP', 'NC_010473.1')}
		
	# ---------------------------------- #
	
	reference_length = { 
			'RNASeq_Control':(5455, 4576, 4686137),
			'RNASeq_State_1':(11792, 4311, 4686137),
			'RNASeq_State_2':(11792, 4311, 4686137),
			'RNASeq_State_3':(11792, 4311, 4686137),
			'RNASeq_State_4':(11792, 4311, 4686137),
			'RNASeq_State_5':(11792, 4311, 4686137),
			'RNASeq_State_6':(11792, 4311, 4686137),
			'RNASeq_State_7':(11792, 4311, 4686137),
			'RNASeq_State_8':(11792, 4311, 4686137),
			'RiboSeq_Control':(5455, 4576, 4686137),
			'RiboSeq_State_1':(11792, 4311, 4686137),
			'RiboSeq_State_2':(11792, 4311, 4686137),
			'RiboSeq_State_3':(11792, 4311, 4686137),
			'RiboSeq_State_4':(11792, 4311, 4686137),
			'RiboSeq_State_5':(11792, 4311, 4686137),
			'RiboSeq_State_6':(11792, 4311, 4686137),
			'RiboSeq_State_7':(11792, 4311, 4686137),
			'RiboSeq_State_8':(11792, 4311, 4686137)}
		
	# ---------------------------------- #
	
	HiSeq_list = "list of the names of the raw FASTQ files"
	
	# ---------------------------------- #
	
	sample_list = {}
	for HiSeq_name in HiSeq_list:
		sample_list[HiSeq_name] = "sample_name"		# for example --> sample_list['RC'] = 'RNASeq_Control'
	
	# ---------------------------------- #
	
	Linker_Sequence = 'CTGTAGGCACCATCAATATCTCGTATGCCGTCTTCTGCTTG'
	
	# ---------------------------------- #
	
	RAW_FASTQ_filename = {}
	Trimmed_RAW_FASTQ_filename = {}
	Unaligned_FASTQ_filename = {}
	Coding_Genes_Fasta_filename = {}
	Coding_Genes_INDEXED_Fasta_Filename = {}
	Coding_Genes_GFF_filename = {}
	Coding_Genes_MAP_output_filename = {}
	circuit_output_RD_exponential_fit_file = {}
	reporter_output_RD_exponential_fit_file = {}
	genome_output_RD_exponential_fit_file = {}
	
	for HiSeq_name in HiSeq_list:
		RAW_FASTQ_filename[sample_list[HiSeq_name]] = "directory for the input raw FASTQ files" + "/" + HiSeq_name + ".fastq"
		Trimmed_RAW_FASTQ_filename[sample_list[HiSeq_name]] = "directory for the input raw FASTQ files" + "/" + HiSeq_name + "_trimmed.fastq"
		Unaligned_FASTQ_filename[sample_list[HiSeq_name]] = "directory for the output unaligned FASTQ files" + "/" + HiSeq_name + "_unaligned.fastq"
		Coding_Genes_Fasta_filename[sample_list[HiSeq_name]] = "directory for the input FASTA files" + "/" + sample_list[HiSeq_name] + ".fasta"
		Coding_Genes_INDEXED_Fasta_Filename[sample_list[HiSeq_name]] = "directory for the output INDEXED FASTA files" + "/" + sample_list[HiSeq_name]
		Coding_Genes_GFF_filename[sample_list[HiSeq_name]] = "directory for the input GFF files" + "/" + sample_list[HiSeq_name] + ".gff"
		Coding_Genes_MAP_output_filename[sample_list[HiSeq_name]] = "directory for the output map files" + "/" + sample_list[HiSeq_name] + ".map"
		circuit_output_RD_exponential_fit_file[sample_list[HiSeq_name]] = "directory for the RD correction fitting plots" + "/" + sample_list[HiSeq_name] + "_circuit.pdf"
		reporter_output_RD_exponential_fit_file[sample_list[HiSeq_name]] =  "directory for the RD correction fitting plots" + "/" + sample_list[HiSeq_name] + "_reporter.pdf"
		genome_output_RD_exponential_fit_file[sample_list[HiSeq_name]] =  "directory for the RD correction fitting plots" + "/" + sample_list[HiSeq_name] + "_genome.pdf"
		
	# ---------------------------------- #
	
	tRNA_rRNA_operon_GFF_filename = "directory for the GFF file of tRNA-rRNA operon positions"
	
	output_RNASeq_reads_filename = "directory for the output mapped RNASeq reads files"
	output_RiboSeq_reads_filename = "directory for the output mapped RiboSeq reads files"
	output_RNA_seq_Profiles_filename = "directory for the output RNASeq profile file"
	output_Ribo_seq_Profiles_filename = "directory for the output RiboSeq profile file"
	
	output_FPKM_table_filename = "directory for the output FPKM table file"
	output_RD_table_filename = "directory for the output RD table file"
	output_Protein_MW_table_filename = "directory for the output Protein MW table file"
	output_Proteome_Fraction_table_filename = "directory for the output Proteome Fraction table file"
	
	# ---------------------------------- #
	
	test = Genetic_Circuit_Analysis_using_Sequening_Methods(par, inducer, degradation, RNAP_p15A_total, RNAP_p15A, coef1, coef2, omega, total_active_ribosomes, codon_table, amino_acid_MW, frag_dist_bin, \
								reference_name, reference_length, Linker_Sequence, RAW_FASTQ_filename, Trimmed_RAW_FASTQ_filename, Unaligned_FASTQ_filename, \
								Coding_Genes_Fasta_filename, Coding_Genes_INDEXED_Fasta_Filename, Coding_Genes_GFF_filename, Coding_Genes_MAP_output_filename, \
								circuit_output_RD_exponential_fit_file, reporter_output_RD_exponential_fit_file, genome_output_RD_exponential_fit_file, tRNA_rRNA_operon_GFF_filename, \
								output_RNASeq_reads_filename, output_RiboSeq_reads_filename, output_RNA_seq_Profiles_filename, output_Ribo_seq_Profiles_filename, \
								output_FPKM_table_filename, output_RD_table_filename, output_Protein_MW_table_filename, output_Proteome_Fraction_table_filename)
	
	# ---------------------------------- #
	
	test.map_raw_FASTQ_files()
	
	# ---------------------------------- #
	
	RNA_seq_Read, Ribo_seq_Read = test.generate_mapped_reads()
	
	# ---------------------------------- #
	
	RNA_seq_Profile, RNA_seq_total_read, RNA_seq_total_read_tRNA_rRNA_operon_deleted, RNA_seq_total_tRNA_fraction = test.generate_RNA_Seq_profile()
	
	# ---------------------------------- #
	
	Ribo_seq_Profile, Ribo_seq_total_read, Ribo_seq_total_read_tRNA_rRNA_operon_deleted, Ribo_seq_total_tRNA_fraction = test.generate_Ribo_Seq_profile()
	
	# ---------------------------------- #
	
	fragment_size_distribution = test.generate_mapped_fragment_size_distribution()
	
	# ---------------------------------- #
	
	sample_name = 'RNASeq_State_1'		# this is an example
	Chromosome = '0x41v70'
	Borders = [1910, 8780]
	alfa  = 2.25		# 2.25 for circuit plasmid
	#alfa  = 1.0		# 1.0 for reporter plasmid
	diff_log_pos, diff_log_neg_rev, TSS_fwd, TSS_rev = test.identify_promoter_TSS(RNA_seq_Profile, sample_name, Chromosome, Borders, alfa)
	
	# ---------------------------------- #
	
	sample_name = 'RNASeq_State_1'		# this is an example
	Chromosome = '0x41v70'
	Strand = 'fwd'
	X2 = 20
	X1 = 10
	alfa  = 2.25		# 2.25 for circuit plasmid
	#alfa  = 1.0		# 1.0 for reporter plasmid
	promoter_activity = {}
	for TSS in TSS_fwd[sample_name][Chromosome]:
		promoter_activity[TSS] = test.calculate_promoter_activity(RNA_seq_Profile, sample_name, Chromosome, Strand, TSS, X2, X1, alfa)
	
	# ---------------------------------- #
	
	Promoter_list_cir = ['pSrpR','pBetI','pBAD1','pTet1','pPhlF','pHlyIIR','pAmtR','pAmeR','pTet2','pTac','pBAD2']
	Promoter_list_rep = ['pBM3R1']
	alfa_cir  = 9.0			# 9 for circuit plasmid
	alfa_rep = 4.0			# 4 for reporter plasmid
	Cello_flux_cir, Cello_flux_rep = test.calculate_Cello_promoter_activity(Promoter_list_cir, Promoter_list_rep, alfa_cir, alfa_rep)	# unit of flux is RNAP/s
	
	# ---------------------------------- #
	
	sample_name = 'RNASeq_State_1'		# this is an example
	Chromosome = '0x41v70'
	Strand = 'fwd'
	ribozyme_cut_position = 2069		# PhlF ribozyme
	offset = 160
	ribozyme_cleavege_efficiency = test.calculate_ribozyme_cleavege_efficiency(RNA_seq_Read, sample_name, Chromosome, Strand, ribozyme_cut_position, offset)
	
	# ---------------------------------- #
	
	sample_name = 'RNASeq_State_1'		# this is an example
	Chromosome = '0x41v70'
	Borders = [1910, 8780]
	alfa  = 2.25		# 2.25 for circuit plasmid
	#alfa  = 1.0		# 1.0 for reporter plasmid
	window_size = 10
	diff_log_pos, diff_log_neg_rev, TTS_window_fwd, TTS_window_rev, second_diff_log, second_diff_log_rev, TTS_fwd, TTS_rev = test.identify_terminator_TTS(RNA_seq_Profile, sample_name, Chromosome, Borders, alfa, window_size)
	
	# ---------------------------------- #
	
	sample_name = 'RNASeq_State_1'		# this is an example
	Chromosome = '0x41v70'
	Strand = 'fwd'
	X2 = 20
	X1 = 10
	alfa  = 2.25		# 2.25 for circuit plasmid
	#alfa  = 1.0		# 1.0 for reporter plasmid
	terminator_strength = {}
	for TTS in TTS_fwd[sample_name][Chromosome]:
		terminator_strength[TTS] = test.calculate_terminator_strength(RNA_seq_Profile, sample_name, Chromosome, Strand, TTS, X2, X1, alfa)
	
	# ---------------------------------- #
	
	sample_name = 'RNASeq_State_1'		# this is an example
	Chromosome = '0x41v70'
	gene = 'PhlF'
	GFF_filename = Coding_Genes_GFF_filename['RNASeq_State_1']
	FPKM = test.calculate_FPKM(GFF_filename, RNA_seq_Profile, sample_name, Chromosome, gene)
	
	# ---------------------------------- #
	
	GFF_filename = Coding_Genes_GFF_filename['RNASeq_State_1']
	FPKM_dict = test.generate_FPKM_table(GFF_filename, RNA_seq_Profile)
	
	# ---------------------------------- #
	
	sample_name = 'RNASeq_State_1'		# this is an example
	Chromosome = '0x41v70'
	gene = 'PhlF'
	alfa  = 2.25		# 2.25 for circuit plasmid
	#alfa  = 1.0		# 1.0 for reporter plasmid
	GFF_filename = Coding_Genes_GFF_filename['RNASeq_State_1']
	mRNA_level = test.calculate_mRNA_level(GFF_filename, RNA_seq_Profile, sample_name, Chromosome, gene, alfa)
	
	# ---------------------------------- #
	
	GFF_filename = Coding_Genes_GFF_filename['RNASeq_State_1']
	aa_ends_excluded = 5
	RD_dict = test.generate_RD_table(GFF_filename, Ribo_seq_Profile, aa_ends_excluded)
	
	# ---------------------------------- #
	
	sample_name = 'RNASeq_State_1'		# this is an example
	Chromosome = '0x41v70'
	gene = 'PhlF'
	alfa  = 2.25		# 2.25 for circuit plasmid
	#alfa  = 1.0		# 1.0 for reporter plasmid
	GFF_filename = Coding_Genes_GFF_filename['RNASeq_State_1']
	TE = test.calculate_translation_efficiency(GFF_filename, RNA_seq_Profile, sample_name, Chromosome, gene, alfa, RD_dict)
	
	# ---------------------------------- #
	
	FASTA_filename = Coding_Genes_Fasta_filename['RNASeq_State_1']
	GFF_filename = Coding_Genes_GFF_filename['RNASeq_State_1']
	Protein_MW = test.calculate_MW(FASTA_filename, GFF_filename)
	
	# ---------------------------------- #
	
	test.generate_MW_table(Protein_MW)
	
	# ---------------------------------- #
	
	Proteome_Fraction = test.calculate_proteome_fraction(RD_dict, Protein_MW)
	
	# ---------------------------------- #
	
	test.generate_proteome_fraction_table(Proteome_Fraction)
	
	# ---------------------------------- #
	
	promoter_TSS_cir = {'pSrpR':  1976,'pBetI':  2070,'pBAD1':  3126,'pTet1':  3236,'pPhlF':  4142,'pHlyIIR':4218,'pAmtR':  5027,'pAmeR':  5103,'pTet2':  5913,'pTac':   6809,'pBAD2':  7922}
	promoter_TSS_rep = {'pBM3R1':1030}
	Gates = { 'PhlF':['pSrpR','pBetI'], 'SrpR':['pBAD1','pTet1'],'BM3R1':['pPhlF','pHlyIIR'],'BetI':['pAmtR','pAmeR'],'AmeR':['pTet2'],'HlyIIR':['pTac'],'AmtR':['pBAD2'],'YFP':['pBM3R1']}
	gate_list = ['PhlF','SrpR','BM3R1','BetI','AmeR','HlyIIR','AmtR']
	RBSs = [(2156,2156+150),(3342,3342+150),(4326,4326+150),(5189,5189+150),(6016,6016+150),(6886,6886+150),(8024,8024+150),(1124,1124+150)]
	terminators = [(2759,2759+57),(3984,3984+90),(4905,4905+57),(5777,5777+47),(6676,6676+47),(7522,7522+53),(8693,8693+57),(1844,1844+61)]
	alfa_1 = 2.25		# 2.25 for circuit plasmid
	alfa_2 = 1.0		# 1.0 for reporter plasmid
	X2 = 20
	X1 = 10
	Gate_input_Cello, Gate_output_Cello, Gate_input_RNAP_Flux, Gate_output_RNAP_Flux, Gate_ribosome_usage = test.generate_gate_response_function(RNA_seq_Profile, Ribo_seq_Profile, alfa_1, alfa_2, promoter_TSS_cir, promoter_TSS_rep, X2, X1, Gates, RBSs, terminators, gate_list)
	
	# ---------------------------------- #
	
	alfa_1 = 2.25		# 2.25 for circuit plasmid
	alfa_2 = 1.0		# 1.0 for reporter plasmid
	border_circuit_control = (1910,5316)
	border_reporter_control = (892,2169)
	border_circuit = (1910,11653)
	border_reporter = (892,1904)
	Total_RNAP_Flux = test.calculate_circuit_total_RNAP_flux(RNA_seq_Profile, alfa_1, alfa_2, border_circuit_control, border_reporter_control, border_circuit, border_reporter)
	
	# ---------------------------------- #
	
	border_circuit_control = (1910,5316)
	border_reporter_control = (892,2169)
	border_circuit = (1910,11653)
	border_reporter = (892,1904)
	Total_Ribo_usage = test.calculate_total_ribosome_usage(Ribo_seq_Profile, border_circuit_control, border_reporter_control, border_circuit, border_reporter)
	
	# ---------------------------------- #
	
