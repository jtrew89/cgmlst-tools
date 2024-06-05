#!/usr/bin/python

##Import modules
import pandas as pd
import argparse

def main(args):
	##load data
	alle_prof_df = pd.read_table(args['pro_db']).set_index('FILE')
	
	##Variables used in script
	isolates = list(alle_prof_df.index)
	out_alle_df = pd.DataFrame(index=isolates, columns=isolates)
	
	##Replace acronyms from chewbac output and treat as missing data '0'
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
	##loop throgh getting pairwise distance (allele difference) for all isolates
	for isolate_1 in isolates:
	    
	    for isolate_2 in isolates:
	        dist = len(alle_prof_df.loc[isolate_1].compare(alle_prof_df.loc[isolate_2]))
	        out_alle_df.at[isolate_1, isolate_2] = dist
	
	
	final_alle_df.to_csv(args['out_dir'], sep='\t', index=True)
 
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
	args = parser.parse_args()

	##Convert the argparse.Namespace to a dictionary: vars(args)
	main(vars(args))
	sys.exit(0)
