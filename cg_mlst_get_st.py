#!/usr/bin/env python

##Load modules
import pandas as pd
import argparse
import time
import re
import sys
pd.options.mode.chained_assignment = None  # default='warn'

##script in function so that it can be used as a module in the main run_snippy.py
def main(args):

	##Load dfs used in script
	if str(args['in_query']).endswith('.csv'):
		query_df = pd.read_csv(args['in_query'], low_memory=False)
	else:
		query_df = pd.read_table(args['in_query'], low_memory=False)
	query_df.set_index('FILE', inplace=True) ##Set isolate id as index so that loopin can be used on index
	query_df = query_df.applymap(str)
	if str(args['db']).endswith('.csv'):
		profiles_df = pd.read_csv(args['db'], low_memory=False)
	else:
		profiles_df = pd.read_table(args['db'], low_memory=False)

	profiles_df.set_index('ST', inplace=True)
	profiles_df = profiles_df.applymap(str) #make df str otherwise there is issue with replace

	##Variables used in script
	isolate_id = list(query_df.index)
	sts_results = []
	isolate_results = []
	novel_sts = []
	novel_isolates = []
	novel_st_multi = []
	novel_iso_multi = []

	##Replace acronyms from chewbac output and treat as missing data '0'
	query_df.replace(
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

	##Loop though query (using isolate ID), and check profiles against profile reference to get ST
	for isolate in isolate_id:

		current_st = profiles_df.iloc[-1].name #get current latest st
		current_st = re.sub('_','', str(current_st)) #remove novel st identifier (if present)
		novel_st = str(int(current_st) + 1) #next novel st
		novel_profile = query_df.loc[isolate] #current allele profile if novel
		query_results = profiles_df[profiles_df == query_df.loc[isolate]].dropna() #conditional selection for rows that match exactly (get ST)

		if not query_results.empty: #if the datafame is not empty, that means there is an exact match and it is not a novel ST in this run
			st = query_results.iloc[0].name #get the matchin ST
			st_nov_check = re.sub('_','', str(st)) #variable to check novelty of st
			if str(st_nov_check) in novel_sts: #check to see if any of the sts in sts_results were novel during this run
				novel_st_multi.append(st)
				novel_iso_multi.append(isolate)
			else:
				sts_results.append(st) #append ST to list
				isolate_results.append(isolate) #append isolate seached in list
		elif query_results.empty: ##if the datafame is empty, that means there is not a match and the profile is novel (new ST)
			novel_sts.append(novel_st) #appened novel ST to list
			novel_isolates.append(isolate) #appened isolate with novel allele profile to list
			#novel_profile['ST'] = '_'+novel_st #add novel ST to novel profile series (perhaps no longer needed. Review later)
			novel_profile.rename('_'+novel_st, inplace=True) #rename novel series with novel ST
			profiles_df = profiles_df.append(novel_profile) #append novel profile to reference profile db

	##Make table of isolates and their STs (nevel or not)
	run_novel_sts_out = pd.DataFrame({'Isolate_ID':novel_isolates,'ST':novel_sts})
	run_novel_sts_out.to_csv(args['out_dir']+'run_novel_sts_out.csv',index=False)

	##Variables for stats
	tot_isos = len(isolate_id)
	tot_nov = len(novel_sts)
	tot_iden = len(sts_results) - len(novel_st_multi)
	tot_re_st = len(novel_st_multi)

	if not sts_results:
		print('All isolates had novel STs')
	else:
		run_ident_sts_out = pd.DataFrame({'Isolate_ID':isolate_results,'ST':sts_results})
		run_ident_sts_out.to_csv(args['out_dir']+'run_ident_sts_out.csv',index=False)
		print(f'Of {tot_isos} isolate(s), {tot_nov} had novel ST(s), {tot_iden} had existing ST(s) and {tot_re_st} isolate(s) were identical in this run')

	##Save new reference db
	timestr = time.strftime("%d_%m_%Y-%H_%M")
	profiles_df.to_csv(f"{args['db_dir']}profiles.list.{timestr}.csv")

if __name__ == '__main__':
	##Create arguments
	parser = argparse.ArgumentParser(description='Script written to get STs for isolates just run through chebbaca and output/add to reference db of entero allele schema')
	parser.add_argument('-i', '--input', dest='in_query', help='Name of and path to input file, (Has to be comma seperated)', required=True)
	parser.add_argument('-db', '--profiles_db', dest='db', help='Path to profiles.list (whichever version) and file name')
	parser.add_argument('-dd', '--db_dir', dest='db_dir', help='Directory to output update db', required=True)
	parser.add_argument('-od', '--output_dir', dest='out_dir', help='Directory for output file containing list of isolates run and the assigned STs (novel to this run:run_novel_sts_out.csv, or alreadr present in the version of the schema database used:run_novel_sts_out.csv)', required=True)

	args = parser.parse_args()

	##Convert the argparse.Namespace to a dictionary: vars(args)
	main(vars(args))
	sys.exit(0)
