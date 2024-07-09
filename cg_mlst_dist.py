#!/usr/bin/env python

##Import modules
import pandas as pd
import argparse
import numpy as np
import sys

def main(args):
	##load data
	alle_prof_df = pd.read_table(args['pro_db']).set_index('Isolate_ID')

	##Variables used in script
	isolates = list(alle_prof_df.index)
	out_alle_df = pd.DataFrame(index=isolates, columns=isolates)

	##Replace acronyms from chewbac and string from etoki output and treat as missing data '0'
	loci = alle_prof_df.columns[1:]
	for locus in loci:
		alle_prof_df[locus] = (pd.to_numeric(alle_prof_df[locus], errors='coerce'). fillna(0))

	'''
	alle_prof_df.replace(
	   		[
	   		'INF-','LNF',
	   		'PLNF', 'PLOT3',
	   		'PLOT5', 'LOTSC',
	   		'NIPHEM', 'NIPH',
	   		'PAMA', 'ALM', 'ASM'
	   		],
	   		[
	   		'', '0','0','0','0','0',
	   		'0','0','0','0','0'
	   		],
	   		inplace=True, regex=True
	   			)
	'''

	##loop throgh getting pairwise distance (allele difference) for all isolates
	if args['miss']:
		for isolate_1 in isolates:

			for isolate_2 in isolates:
				dist = len(alle_prof_df.loc[isolate_1].compare(alle_prof_df.loc[isolate_2]).replace(0.0,np.nan).dropna()) #drops missing data per-pairwise comparison
				out_alle_df.at[isolate_1, isolate_2] = dist


		out_alle_df.to_csv(args['out_dir'] + 'allele_dist.tsv', sep='\t', index=True)
	else:
		for isolate_1 in isolates:

			for isolate_2 in isolates:
				dist = len(alle_prof_df.loc[isolate_1].compare(alle_prof_df.loc[isolate_2])) #keeps missing data
				out_alle_df.at[isolate_1, isolate_2] = dist

		out_alle_df.to_csv(args['out_dir'] + 'allele_dist.tsv', sep='\t', index=True)

if __name__=='__main__':
	##Create arguments
	parser = argparse.ArgumentParser(description='Script written to create distance matrix from results_alleles.tsv, of some other tab-delimeted allele profile')

	parser.add_argument(
	'-pdb', '--profile_db', dest='pro_db',
	help='Alelle profile db with isolates you wish analysed (Include path)',
	required=True
	)
	parser.add_argument(
	'-od', '--output_dir', dest='out_dir',
	help='Directory to output matrix',
	required=True
	)
	parser.add_argument(
	'-m', '--missing_data', dest='miss',
	help='Missing data ("0"s) will not be included in pairwise comparisons',
	required=False, action='store_true'
	)
	args = parser.parse_args()

	##Convert the argparse.Namespace to a dictionary: vars(args)
	main(vars(args))
	sys.exit(0)
