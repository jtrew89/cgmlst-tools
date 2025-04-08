#!/usr/bin/env python

##Import modules
import pandas as pd
import argparse
import numpy as np
import sys
import numpy as np
from tqdm import tqdm


def main(args):
	##load data
	alle_prof_df = pd.read_table(args['pro_db']).set_index('Isolate_ID')
	#alle_prof_df = alle_prof_df[~alle_prof_df.index.duplicated(keep='first')] ##remove duplicates

	##Variables used in script
	isolates = list(alle_prof_df.index)
	loci_num = len(alle_prof_df.columns)

	##Replace acronyms from chewbac and string from etoki output and treat as missing data '0'
	loci = alle_prof_df.columns[1:]
	for locus in loci:
		alle_prof_df[locus] = (pd.to_numeric(alle_prof_df[locus], errors='coerce'). fillna(0))

	"""
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
	"""

	##remove missing data based on thresholds
	if args['miss']:
		alle_prof_df.replace(0, np.nan,inplace=True) #to act upon missing data
		isolates_num = len(alle_prof_df.index)
		alle_prof_df.dropna(thresh=len(alle_prof_df.columns) - len(alle_prof_df.columns)/100 * int(args['i_threshold']),inplace =True) #drop isolates with missing data exceeding threshold
		alle_prof_df = alle_prof_df.dropna(axis='columns') #remove loci with missing data from dataframe
		miss = loci_num - len(alle_prof_df.columns) #number of droped loci due to missing data
		droped_isolates = isolates_num - len(alle_prof_df.index) #number of dropped isolates due to missing data
		isolates = list(alle_prof_df.index)
		out_alle_df = pd.DataFrame(index=isolates, columns=isolates) #empy df for output
		if miss >= loci_num/100 * int(args['m_threshold']): #if number of dropped loci is greater than threshold set, exit script
			sys.exit(f"Loci with missing data exceeds missing data threshold ({miss} of {loci_num} loci with missing data and {droped_isolates} of {isolates_num} isolates dropped, due to individual % of missing data)")

		##loop throgh getting pairwise distance (allele difference) for all isolates
		for isolate_1 in tqdm(isolates):

			for isolate_2 in isolates:
				dist = len(alle_prof_df.loc[isolate_1].compare(alle_prof_df.loc[isolate_2]))
				out_alle_df.at[isolate_1, isolate_2] = dist
			out_alle_df.to_csv(args['out_dir'] + 'allele_dist_miss.tsv', sep='\t', index=True)
		print(str(miss) + " loci ignored, due to missing '0s' data exceeding threshold")
		print(str(droped_isolates) + " isolates removed, due to missing '0s' data exceeding threshold")

	else:
		out_alle_df = pd.DataFrame(index=isolates, columns=isolates)
		for isolate_1 in tqdm(isolates):

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
	parser.add_argument(
	'-mt', '--missing_threshold', dest='m_threshold',
	help='Percenatage of total missing data accepted, across all isolates (putting 5 will mean no less than 95 %% of data present). If lower, the process will stop',
	required=False, default=10, type=int
	)
	parser.add_argument(
	'-it', '--isolate_threshold', dest='i_threshold',
	help='Percentage of missing data accepted per isolate (putting 5 will mean no less than 95 %% of data present). Isolates not meeting this percentage will be removed',
	required=False, default=10, type=int
	)

	args = parser.parse_args()

	##Convert the argparse.Namespace to a dictionary: vars(args)
	main(vars(args))
	sys.exit(0)
