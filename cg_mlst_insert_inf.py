#!/usr/bin/env python

"""
Script written to take the novel_alleles.fasta output file from chewBBACA and
insert each novel allele into their respective schema with a 'inf_' mark. So that
novel alleles are maked int he schema, as chewBBACA currently doesnt mark novel
alleles
"""

##Import modules
import os
import argparse
import pandas as pd
import sys
import re

##Make whole script into function so that the script can be imported by another script
def main(args):

	#Variables for input
	novel_fas = args['in_dir']+'novel_alleles.fasta'  # Input fasta file
	sch_fas = args['schm_fas'] # Fatsta schema used by etoki

	with open(novel_fas, 'r') as file:
		novel_lines = file.readlines()

	#Function to find latest allele for a give loci, so the novel sequences can be inseryter beneath
	def novel_alle_num(loci):
		schema = open(sch_fas)
		schema_read = schema.readlines()
		alle_n_list = []
		for line in schema_read:
			if loci in line:
				alle_n_list.append(int(line.split('_')[2]))
		return max(alle_n_list)

	#Loop through each novel sequence and insert the sequence and header into schema fasta file
	for num, line in enumerate(novel_lines):
		if '>' in line:
			curr_loci = re.sub('>','',line.rsplit('_', 2)[0]) #remove unique allele number to use to send to the out fasta file
			header_insert = line.rsplit('_', 1)[0] + '_' + '\n'
			header_search = curr_loci + '_' + str(novel_alle_num(curr_loci)) #
			alle_seq = novel_lines[num + 1]
			seq_insert =  alle_seq
			with open(sch_fas, 'r') as file:
				sch_lines = file.readlines()
			header_insert_line = ''
			seq_insert_line = ''
			for num, line in enumerate(sch_lines):
				if header_search in line:
					header_insert_line = num + 2
					seq_insert_line = num +3
				else:
					pass
			#test to see is the search function is working, used in testing
			if header_insert_line == '':
				sys.exit(f'{header_search} allele not found in search')

			sch_lines.insert(header_insert_line, header_insert)
			sch_lines.insert(seq_insert_line, seq_insert)

			with open(args['schm_fas'], 'w') as file:
				file.writelines(sch_lines)

		else:
			pass
if __name__ == '__main__':
	##Create arguments
	parser = argparse.ArgumentParser(description='Script written to take the novel_alleles.fasta output file from chewBBACA and insert each novel allele into their respective schema with a "inf_" mark. So that novel alleles are maked in the schema, as chewBBACA currently doesnt mark novel alleles.')
	parser.add_argument(
			'-id', '--input_directory',
			dest='in_dir',
			help='Directory where novel_alleles.fasta is kept',
			required=True
			)
	parser.add_argument(
			'-s', '--scheema_fasta',
			dest='schm_fas',
			help='fasta file with 3002 schema',
			required=True
			)

	args = parser.parse_args()

	##Convert the argparse.Namespace to a dictionary: vars(args)
	main(vars(args))
	sys.exit(0)
