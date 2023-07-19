#!/usr/bin/env python

"""
Script written to take the novel_alleles.fasta output file from chewBBACA and
insert each novel allele into their respective schema with a 'inf_' mark. So that
novel alleles are maked int he schema, as chewBBACA currently doesnt mark novel
alleles
"""

##Import modules
from Bio import SeqIO
import os
import argparse
import pandas as pd
import sys

##Make whole script into function so that the script can be imported by another script
def main(args):

	##Some code found online that I should be able to change into what I want
	fasta_file = args['in_dir']+'novel_alleles.fasta'  # Input fasta file
	id_file =  args['in_dir']+'novel_alleles.csv' # Input novel allele ID
	loci_dir = args['schm_dir'] # Output fasta file which allele is sent to

	id_df = pd.read_csv(id_file)
	id_df['fas_id'] = id_df.Locus + '_' + id_df['Novel Allele'].apply(str)
	wanted = sorted(set(id_df['fas_id']))

	#fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')  #if I want to iterate through fasta sequences at a later date
	fasta_dict = SeqIO.index(args['in_dir']+'novel_alleles.fasta', 'fasta') #Create a dictionary out of the fasta so that is can be indexed for a particular sequence

	end = False
	for alle in wanted:
		fas_out = alle.rsplit('_', 1)[0] #remove unique allele number to use to send to the out fasta file
		with open(loci_dir+fas_out+'.fasta', "a") as f:
			curr_seq = fasta_dict[alle]
			curr_seq.id = fas_out + '_inf_' + curr_seq.id[curr_seq.id.rindex('_')+1:]
			SeqIO.write(curr_seq, f, "fasta")

if __name__ == '__main__':
	##Create arguments
	parser = argparse.ArgumentParser(description='Script written to take the novel_alleles.fasta output file from chewBBACA and insert each novel allele into their respective schema with a "inf_" mark. So that novel alleles are maked in the schema, as chewBBACA currently doesnt mark novel alleles.')
	parser.add_argument(
			'-id', '--input_directory',
			dest='in_dir',
			help='Directory where novel_alleles.fasta and novel alleles.csv are kept',
			required=True
			)

	parser.add_argument(
			'-sd', '--schema_directory',
			dest='schm_dir',
			help='Directory where 3002 loci fasta files are kept',
			required=True
			)
	args = parser.parse_args()

	##Convert the argparse.Namespace to a dictionary: vars(args)
	main(vars(args))
	sys.exit(0)
