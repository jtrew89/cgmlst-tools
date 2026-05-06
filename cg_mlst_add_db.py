#!/usr/bin/env python

##Load modules
import pandas as pd
import os
import argparse
import time
import re
import sys

def main(args):
	##Variables used in script
	date = time.strftime('%d_%m_%Y')
	in_file = args['in_dir']+'results_alleles.tsv'
	db_file = os.path.basename(args['str_db'])
	db_dir = re.sub(db_file,'',args['str_db']) #get db dir from path to db by removing input db name from path

	##Load database and allele profile to be added
	in_alle = pd.read_table(in_file, low_memory=False)
	db_alle = pd.read_csv(args['str_db'], low_memory=False)

	##Concat allele profiles for new db
	new_db = pd.concat([in_alle,db_alle], ignore_index=True)

	##Save new db
	new_db.to_csv(db_dir+'curr_db_'+date+'.tsv', sep='\t', index=False)

if __name__ == '__main__':
	##Create arguments
	parser = argparse.ArgumentParser(description="Script written to add a new results_alleles.tsv allele profile to the Salmonella enterica database (sent_curr_db.csv)")
	parser.add_argument(
		'-id', '--input_directory',
		dest='in_dir',
		help="Directory where results_alleles.tsv is kept",
		required=True
			)
	parser.add_argument(
		'-sdb','--strain_db',
		dest='str_db',
		help="Strain db to update and path to",
		required=True
			)

	args = parser.parse_args()

	##Convert the argparse.Namespace to a dictionary: vars(args)
	main(vars(args))
	sys.exit(0)