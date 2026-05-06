#!/usr/bin/env python

"""
Script written to identify INF alleles from results_alleles.tsv output from
Chewbbaca and to mark all novel allele designations with...
"""
##Import libraries
import pandas as pd
import argparse
import re
import sys

##Make whole script into function so that the script can be imported by another script
def main(args):
	##Load allele profiles
	al_prof_df = pd.read_table(f"{args['in_dir']}results_alleles.tsv")

	##Creat Variables used in script
	loci = al_prof_df.columns.tolist()
	locus_alle = []
	novel_alle = []
	novel_df = pd.DataFrame()

	##Loop through df and identify cells with INF
	for locus in loci:
		curr_loci = al_prof_df.loc[:,locus]

		for alle in curr_loci:
			if "_" in str(alle):
				alle = re.sub('_','',alle)
				novel_alle.append(alle)
				locus_alle.append(locus)
			else:
				pass

	##Create df to group Loci and nevel allele information

	novel_df['Locus'] = locus_alle
	novel_df['Novel Allele'] = novel_alle

	##Output novel alleles
	novel_df = novel_df.drop_duplicates(subset = ['Locus','Novel Allele'])
	novel_df.to_csv(f"{args['out_dir']}novel_alleles.csv", index=False)

	## Old section required for chewbbaca output
	"""
	INF-'novel allele' is only statted once, any other isolates that
	have that allele in that run will just have the allele number
	(e.g INF-325 for the first isolates with the novel allele,
	325 for every isolate after)


	##Add mark to all novel alleles in results_alleles.tsv
	form_results = pd.DataFrame()

	for locus in loci:
		curr_loci = al_prof_df.loc[:,locus]

		for i in curr_loci:
			if "INF" in str(i):
				inf_alle = re.sub('INF-', '', i)
				curr_loci = curr_loci.str.replace('INF-'+inf_alle, inf_alle)
				curr_loci = curr_loci.str.replace(inf_alle, 'inf_'+inf_alle)
				form_results[curr_loci.name] =   curr_loci
		else:
			form_results[curr_loci.name] =   curr_loci

	##Output results_alleles with all novel alleles marked
	form_results.to_csv(f"{args['out_dir']}results_alleles_form.csv", index=False)
	"""

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
