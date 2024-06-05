#!/usr/bin/python

##Import modules
import pandas as pd
import argparse

def main(args):
	##load data
	alle_prof_df = pd.read_table('C:/Users/m1011068/OneDrive - Defra/Documents/Publications/Infantis for Shaun/New_report_tree/all_results_alleles_form_new.tsv').set_index('FILE')
	#alle_prof_df = pd.read_csv('C:/Users/m1011068/OneDrive - Defra/Documents/Publications/Infantis for Shaun/New_report_tree/test_profile.csv').set_index('FILE')
	
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
	
	
	#final_alle_df.to_csv('C:/Users/m1011068/OneDrive - Defra/Documents/Publications/Infantis for Shaun/new_report_tree/distance_mat.tsv',sep='\t')
 
if __name__=='__main__':
	##Create arguments
	parser = argparse.ArgumentParser(description='Script written to identify INF alleles from results_alleles.tsv output from Chewbbaca and to mark all novel allele designations within results_alleles.tsv')
	parser.add_argument(
		'-id', '--input_directory',
		dest='in_dir',
		help="Directory where results_alleles.tsv is kept",
		required=True
                )

	parser.add_argument(
		'-od', '--output_directory',
		dest='out_dir',
		help="""Directory you want the output INF
		list and the formated results_alleles.tsv
		(results_alleles_form.tsv) to be output""",
		required=True
		)
	args = parser.parse_args()

	##Convert the argparse.Namespace to a dictionary: vars(args)
	main(vars(args))
	sys.exit(0)
