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
				for line_2 in con_table_read:
					if line_1.split(',')[1] in line_2:
						alle_n_list.append(int(re.sub('\n','',line_2.split(',')[2])))
				return True, line_1.split(',')[1],line_1.split(',')[2], max(alle_n_list)
			else:
				return False, False, False, False

	##Variables used in script
	os.chdir(args['in_dir'])
	isolates_dic = {}

	##Loop through list of etoki output and pull out allele profile
	novel_alle = {}
	for isolate in glob.glob('*results_alleles.fasta'):
		profile_dic ={}
		isolate_id = re.sub('_results_alleles.fasta','',isolate)
		with open(isolate) as isolate_etoki_out:
			for line in isolate_etoki_out:
				##Checking to see if etoki has given an allele num or an mdtsum (which means the allele is novel or duplicate)
				if '>' in line:
					alle_num = re.sub('id=','',line.split(' ')[2])
					if '-' in alle_num: #if '-' is in the, it is an md5sum
						if md5_finder(alle_num)[0]:
							profile_dic[line.split(' ')[0].split('>')[1]] = int(re.sub('\n','',md5_finder(alle_num)[2]))
						#else:
						#	profile_dic[line.split(' ')[0].split('>')[1]] = '_' + (int(md5_finder(alle_num)[3] + 1))
					else:
						profile_dic[line.split(' ')[0].split('>')[1]] = alle_num
				else:
					pass
		isolates_dic[isolate_id] = profile_dic

	##Put into dataframe
	out_profile = pd.DataFrame(isolates_dic).transpose()
	out_profile.to_csv(args['in_dir'] + 'results_alleles.tsv',sep='\t')

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
