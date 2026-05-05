#!/usr/bin/env python

##Load modules
import pandas as pd
import os
import argparse
import re
import sys
import glob
from tqdm import tqdm

## Main section of script
def main(args):

	## Functions in script
	# Function to check if md5sum from etoki output is in the database (etoki misses some of the duplicates)
	def md5_finder(md5):
		os.chdir(original_dir)
		con_table = open(os.path.abspath(args['con_tab']))
		con_table_read = con_table.readlines()
		for line_1 in con_table_read:
			if md5 == line_1.split(',')[0]:
				return True, line_1.split(',')[1],line_1.split(',')[2]
			elif md5 != line_1.split(',')[0]:
				pass

	# Function to get to ascertain the number for a novel allele
	def novel_alle_num(loci):
		os.chdir(original_dir)
		con_table = open(os.path.abspath(args['con_tab']))
		con_table_read = con_table.readlines()
		alle_n_list = []
		for line in con_table_read:
			if loci in line:
				alle_n_list.append(int(re.sub('_','',re.sub('\n','',line.split(',')[2]))))
		return max(alle_n_list)

	# Function to add novel allele md5 to convert table database
	def convert_table_update(fheader,novel_alle_,loci_d):
		line_search = loci_d + ',' + str(novel_alle_ - 1)
		line_insert = re.sub('value_md5=','',fheader.split(' ')[1]) + ',' + loci_d + ',' + str(novel_alle_) + '_' +'\n'
		line_ind = ''

		os.chdir(original_dir)
		with open(os.path.abspath(args['con_tab']), 'r') as file:
			lines = file.readlines()
		for num, line in enumerate(lines):
			if line_search in line:
			#if line_search + '\n' == line.split(',')[1] + ',' + line.split(',')[2]: #more accurate, but not working right now
				line_ind = num + 1
			else:
				pass

		 #test to see is the search function is working, used in testing
		if line_ind == '':
			sys.exit(f'{line_search} allele not found in search')


		lines.insert(line_ind,line_insert)

		with open(os.path.abspath(args['con_tab']), 'w') as file:
			file.writelines(lines)

	# Function to pick allele if duplicatd using various thresholds
	def start_codon_dup(loci_d,atg,isolate_etoki,md5s,out_count,c_allele_raw):
		current_allele = ''
		#print(loci_d,atg,md5s,out_count,c_allele_raw) # Troubleshooting to test func
		try:
			c_allele = re.sub('_','',c_allele_raw)
		except:
			c_allele = c_allele_raw

		# Set of conditions that decide which of the duplicate alleles are used in the allele profile
		if 'dup_1' in loci_d and isolate_etoki[out_count - 1].startswith('ATG') or 'dup_1' in loci_d and isolate_etoki[out_count - 1].startswith('GTG'):
			atg = int(atg) + 1
			#print('Worked_1') # Troubleshooting 
			if isolate_etoki[out_count + 1].startswith('ATG') or isolate_etoki[out_count + 1].startswith('GTG'):
				atg = int(atg) + 1
				current_allele = min(int(md5s), int(c_allele))
			else:
				current_allele = c_allele
		elif int(atg) > 0 and isolate_etoki[out_count + 1].startswith('ATG') or int(atg) > 0 and isolate_etoki[out_count + 1].startswith('GTG'):
			current_allele = min(int(md5s), int(c_allele))
			#print('Worked_2') # Troubleshooting
		elif int(atg) == 0 and isolate_etoki[out_count + 1].startswith('ATG') or int(atg) == 0 and isolate_etoki[out_count + 1].startswith('GTG'):
			#print('Worked_3') # Troubleshooting
			current_allele = int(md5s)
			atg = int(atg) + 1
		# Not used right now, section for picking allele if neither start a start codon
		elif int(atg) == 0 and not isolate_etoki[out_count + 1].startswith('ATG') or int(atg) == 0 and not isolate_etoki[out_count + 1].startswith('GTG'):
			#print('Worked_4') # Troubleshooting
			current_allele = min(int(md5s), int(c_allele))
		#	if isolate_etoki[out_count + 1].endswith('TAA'):
   		#		current_allele = min(int(md5s), int(c_allele))
		elif int(atg) > 0 and not isolate_etoki[out_count + 1].startswith('ATG') or int(atg) > 0 and not isolate_etoki[out_count + 1].startswith('GTG'):
			current_allele = c_allele
			#print('Worked_5') # Troubleshooting
		#print(atg) # Troubleshooting
		#if int(current_allele) == int(re.sub('_','',c_allele_raw)):
		#	print(c_allele) # tmp test
		#	return c_allele_raw,atg

		#print(current_allele) # Troubleshooting
		return current_allele,atg

 	## Variables used in script
	original_dir = os.getcwd()
	isolates_dic = {}
	novel_alle_out = {}

	## Loop through list of etoki output and pull out allele profile
	for isolate in tqdm(glob.glob(args['in_dir'] + '*results_alleles.fasta')):
		#print(isolate) # Troubleshooting
		profile_dic ={}
		isolate_id = re.sub('_results_alleles.fasta','',isolate)
		isolate_id = re.sub('etoki_out/','',isolate_id)
		isolate_etoki_out = open(isolate)
		isolate_etoki_out_read = isolate_etoki_out.readlines()
		counter = 0 # Check place in convert database
		atg_gtg_start = 0 # For use in start_codon_dup, indicates that previous duplicate of loci has started with atg/gtg, so not to consider any alleles that do not
		for line in isolate_etoki_out_read:
			## Checking to see if etoki has given an allele num or an mdtsum (which means the allele is novel or duplicate)
			if '>' in line:
				alle_num = re.sub('id=','',line.split(' ')[2])
				loci = line.split(' ')[0].split('>')[1]
				if '1d840e66-127e-c54a-1924-4f29a07fb8b8' in alle_num: # This is the md5sum for 'DUPLICATED'. Take duplicated data is treated as missing
					profile_dic[loci] = 0
				elif alle_num.startswith('-'):
					profile_dic[loci] = -1
				elif '-' in alle_num: # If '-' is in the ID, it is an md5sum
					try:
						md5 = str(re.sub('\n','',re.sub('_','',md5_finder(alle_num)[2])))
					except:
						md5 = ''
					loci_dup = re.sub('_dup_.*','',loci) # Remove dup label for downstream parsing if it is a duplicate
					if md5 and 'dup' not in loci: # If allele (md5) is in reference convert_tab and this loci has not been search yet (not a duplicate)
						profile_dic[loci_dup] = md5
						atg_gtg_start = 0
					elif md5 and profile_dic[loci_dup]: # If allele (md5) is in reference convert_tab and this loci has been search before (is a duplicate)
						if args['duplicates']: # If user wants to see duplicate genes for a loci that pass EToKi's thresholds
							profile_dic[loci_dup] = [profile_dic[loci_dup], md5] 
						else:
							profile_dic[loci_dup],atg_gtg_start = start_codon_dup(loci,atg_gtg_start,isolate_etoki_out_read,md5,counter,profile_dic[loci_dup])
					elif not md5 and 'dup' not in loci: # If allele (md5) is not in reference convert_tab (novel allele) and this loci has not been search yet (not a duplicate)
						novel_alle = str(int(novel_alle_num(loci_dup) + 1)) #get novel allele number for current loci
						profile_dic[loci_dup] = '_' + novel_alle
						atg_gtg_start = 0
						if loci_dup in novel_alle_out:
							novel_alle_out[loci_dup].update({novel_alle + '_' + alle_num: isolate_etoki_out_read[counter + 1]})
						else:
							novel_alle_out[loci_dup] = {novel_alle + '_' + alle_num: isolate_etoki_out_read[counter + 1]}
						convert_table_update(line, int(novel_alle),loci_dup) # If md5sum for sequence is not found, allele is novel and the md5sum is addeded to the md5sum database.
					elif not md5 and profile_dic[loci_dup]: # If allele (md5) is not in reference convert_tab (novel allele) and this loci has been search before (is a duplicate)
						novel_alle = str(int(novel_alle_num(loci_dup) + 1)) #get novel allele number for current loci
						if args['duplicates']:
							profile_dic[loci_dup] = [profile_dic[loci_dup], '_' + novel_alle]
						else:
							profile_dic[loci_dup],atg_gtg_start = start_codon_dup(loci,atg_gtg_start,isolate_etoki_out_read,novel_alle,counter,profile_dic[loci_dup])
						if loci_dup in novel_alle_out:
							novel_alle_out[loci_dup].update({novel_alle + '_' + alle_num: isolate_etoki_out_read[counter + 1]})
						else:
							novel_alle_out[loci_dup] = {novel_alle + '_' + alle_num: isolate_etoki_out_read[counter + 1]}
							convert_table_update(line, int(novel_alle),loci_dup)
				else: # If allele has number
					if 'dup' in loci:
						profile_dic[loci] = alle_num
						if isolate_etoki_out_read[counter + 1].startswith('ATG') or isolate_etoki_out_read[counter + 1].startswith('GTG'):
							atg_gtg_start = atg_gtg_start + 1
					elif 'dup' not in loci:
						atg_gtg_start = 0
						profile_dic[loci] = alle_num
			else:
				pass
			counter = counter + 1
		isolates_dic[isolate_id] = profile_dic

		## Put into dataframe
		out_profile = pd.DataFrame(isolates_dic).transpose()
		out_profile.index.name = 'FILE'
		out_profile.to_csv(args['output'], sep='\t')

	##Output novel alleles
	with open(args['in_dir'] + 'novel_alleles.fasta', 'w') as f:
		for key1, value1 in novel_alle_out.items():
			for key2, value2 in value1.items():
				f.write('>' + key1 + '_' + key2 + '\n' + value2)

if __name__=='__main__':
	##Create arguments
	parser = argparse.ArgumentParser(description='A script to convert a file full of etoki fasta outputs "*results_alleles.fasta" in an allele profile for the database')
	parser.add_argument(
		'-id', '--input_directory',
		dest='in_dir',
		help="Directory where *results_alleles.fasta are kept",
		required=True
        )
	parser.add_argument(
		'-ct','--convert_table',
		dest='con_tab',
		help='Convert table used in EtoKi MLST call',
		required=True
		)
	parser.add_argument(
		'-od','--output_file',
		dest='output',
		help='path and name of output file',
		required=True
		)
	parser.add_argument(
		'-d','--show_duplicates',
		dest='duplicates',
		help='Show duplicates called by EToKi in output dataframe',
		required=False,
		action='store_true'
		)
	args = parser.parse_args()

	##Convert the argparse.Namespace to a dictionary: vars(args)
	main(vars(args))
	sys.exit(0)
