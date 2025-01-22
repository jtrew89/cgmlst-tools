#!/usr/bin/env python

##Load modules
import pandas as pd
import os
import argparse
import re
import sys
import glob
import subprocess

##Main section of script
def main(args):

	def md5_finder(md5):

		con_table = open(args['con_tab'])
		con_table_read = con_table.readlines()
		for line_1 in con_table_read:
			alle_n_list = []
			if md5 == line_1.split(',')[0]:
				return True, line_1.split(',')[1],line_1.split(',')[2]
			elif md5 != line_1.split(',')[0]:
				pass

	def novel_alle_num(loci):
		con_table = open(args['con_tab'])
		con_table_read = con_table.readlines()
		alle_n_list = []
		for line in con_table_read:
			if loci in line:
				alle_n_list.append(int(re.sub('\n','',line.split(',')[2])))
		return max(alle_n_list)

	def convert_table_update(fheader,novel_alle_):
		with open(args['con_tab'], 'r') as file:
			lines = file.readlines()

	##Variables used in script
	os.chdir(args['in_dir'])
	isolates_dic = {}

	##Loop through list of etoki output and pull out allele profile
	novel_alle_out = {}

	for isolate in glob.glob('*results_alleles.fasta'):
		profile_dic ={}
		isolate_id = re.sub('_results_alleles.fasta','',isolate)
		isolate_etoki_out = open(isolate)
		isolate_etoki_out_read = isolate_etoki_out.readlines()
		counter = 0
		for line in isolate_etoki_out_read:
			##Checking to see if etoki has given an allele num or an mdtsum (which means the allele is novel or duplicate)
			if '>' in line:
				alle_num = re.sub('id=','',line.split(' ')[2])
				loci = line.split(' ')[0].split('>')[1]
				if '-' in alle_num: #if '-' is in the ID, it is an md5sum
					loci_dup = re.sub('_dup_.*','',loci) #remove dup label for downstream parsing
					novel_alle = str(int(novel_alle_num(loci_dup) + 1)) #get novel allele number for current loci
					if md5_finder(alle_num) and 'dup' not in loci:
						profile_dic[loci_dup] = int(re.sub('\n','',md5_finder(alle_num)[2]))
					elif md5_finder(alle_num) and profile_dic[loci_dup]:
						profile_dic[loci_dup] = [profile_dic[loci_dup], int(re.sub('\n','',md5_finder(alle_num)[2]))]
					elif not md5_finder(alle_num) and 'dup' not in loci:
						profile_dic[loci_dup] = '_' + novel_alle
						novel_alle_out[loci_dup] = {novel_alle + '_' + alle_num: isolate_etoki_out_read[counter + 1]}
						convert_table_update(line, novel_alle)
					elif not md5_finder(alle_num) and profile_dic[loci_dup]:
						profile_dic[loci_dup] = [profile_dic[loci_dup], '_' + novel_alle]
						novel_alle_out[loci_dup] = {novel_alle + '_' + alle_num: isolate_etoki_out_read[counter + 1]}
						convert_table_update(line, novel_alle)
				else:
					profile_dic[loci] = alle_num
			else:
				pass
			counter = counter + 1
		isolates_dic[isolate_id] = profile_dic

	##Put into dataframe
	out_profile = pd.DataFrame(isolates_dic).transpose()
	out_profile.to_csv(args['in_dir'] + 'results_alleles.tsv',sep='\t')

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
	args = parser.parse_args()

	##Convert the argparse.Namespace to a dictionary: vars(args)
	main(vars(args))
	sys.exit(0)
