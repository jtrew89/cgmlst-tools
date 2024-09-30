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

	##Variables used in script
	os.chdir(args['in_dir'])
	isolates_dic = {}

	##Loop through list of etoki output and pull out allele profile
	for isolate in glob.glob('*results_alleles.fasta'):
		profile_dic ={}
		isolate_id = re.sub('_results_alleles.fasta','',isolate)
		with open(isolate) as isolate_etoki_out:
			for line in isolate_etoki_out:
				if '>' in line:
					profile_dic[line.split(' ')[0].split('>')[1]] = re.sub("[^0-9]", "",line.split(' ')[2])
				else:
					pass
		isolates_dic[isolate_id] = profile_dic

	##Put into dataframe
	out_profile = pd.DataFrame(isolates_dic).transpose()

if __name__=='__main__':
	##Create arguments
	parser = argparse.ArgumentParser(description='A script to convert a file full of etoki fasta outputs "*results_alleles.fasta" in an allele profile for the database')
	parser.add_argument(
		'-id', '--input_directory',
		dest='in_dir',
		help="Directory where *results_alleles.fasta are kept",
		required=True
               	)

	args = parser.parse_args()

	##Convert the argparse.Namespace to a dictionary: vars(args)
	main(vars(args))
	sys.exit(0)
