#!/usr/bin/env python

##Load modules
import pandas as pd
import argparse
import time
import re
import sys
import resource
import csv
import os
import shutil
from tqdm import tqdm
pd.options.mode.chained_assignment = None  # default='warn'

# Start timer
start_time = time.time()

##Functions to limit ram
def memory_limit_half():
	"""Limit max memory usage to half."""
	soft, hard = resource.getrlimit(resource.RLIMIT_AS)
	# Convert KiB to bytes, and divide in two to half
	resource.setrlimit(resource.RLIMIT_AS, (get_memory() * 1024 // 2, hard))

def get_memory():
	with open('/proc/meminfo', 'r') as mem:
		free_memory = 0
		for i in mem:
			sline = i.split()
			if str(sline[0]) in ('MemFree:', 'Buffers:', 'Cached:'):
				free_memory += int(sline[1])
	return free_memory  # KiB

##script in function so that it can be used as a module in the main run_snippy.py
def main(args):

	##Chunch size variable (read dataframe in chuncks, as the memory keeps running out)
	chunk_size = 50000

	##Load dfs used in script looped through in chunks
	if str(args['in_query']).endswith('.csv'):
		query_df = pd.read_csv(args['in_query'], low_memory=False)
	else:
		query_df = pd.read_table(args['in_query'], low_memory=False)
	query_df.set_index('FILE', inplace=True) ##Set isolate id as index so that loopin can be used on index
	query_df = query_df.applymap(str)

	##Variables used in script
	isolate_id = list(query_df.index)
	sts_results = []
	isolate_results = []
	novel_profiles_df = pd.DataFrame() #create temporary profile df to store novel profiles
	novel_sts = []
	novel_isolates = []
	novel_st_multi = []
	novel_iso_multi = []
	current_st = ''

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
	for isolate in tqdm(isolate_id):
		novel_profile = query_df.loc[isolate] #current allele profile if novel

		for profiles_chunk in pd.read_csv(args['db'], chunksize=chunk_size, low_memory=False):
			profiles_chunk.set_index('ST', inplace=True)
			profiles_chunk = profiles_chunk.applymap(str) #make df str otherwise there is issue with replace
			novel_multi_results = pd.DataFrame()

			if isolate in isolate_results or isolate in novel_iso_multi: #if there is a result in earlier chunks, don't go through the left over chunks
				pass
			else:
				query_results = profiles_chunk[profiles_chunk == novel_profile].dropna() #conditional selection for rows that match exactly (get ST)
				if not novel_profiles_df.empty:
					if isolate in list(novel_profiles_df.index): #if the current isolate has been added the the novel df already, pass this step
						pass
					else:
						novel_multi_results = novel_profiles_df.drop('ST',axis=1)[novel_profiles_df.drop('ST',axis=1) == novel_profile].dropna() #selecetion to check if the current isolate in the loop matches another isolate from this run, to stop build up of 'novel' STs in the sanme run
				#else:
				#MADE NO SENSE, SO HASHED. WILL LOOK AT IT PROPERLY LATER
				#	novel_multi_results = novel_profiles_df[novel_profiles_df == novel_profile].dropna() #selecetion to check if the current isolate in the loop matches another isolate from this run, to stop build up of 'novel' STs in the sanme run
				if not query_results.empty: #if the datafame is not empty, that means there is an exact match and it is not a novel ST in this run
					st = query_results.iloc[0].name #get the matchin ST
					sts_results.append(st) #append ST to list
					isolate_results.append(isolate) #append isolate seached in list
				elif not novel_multi_results.empty: #if not empty, the profile of this isolate matches that of another in this run and will be asigned the same ST
					nst = re.sub('_','', str(novel_profiles_df.loc[novel_multi_results.index[0]]['ST']))
					novel_st_multi.append(nst)
					novel_iso_multi.append(isolate)
				elif query_results.empty: ##if the datafame is empty, that means there is not a match and the profile is novel (new ST)
					novel_profiles_df = pd.concat([novel_profiles_df, pd.DataFrame([novel_profile])]) #append novel profile to temporary novel profile db
					# DEPRICATED novel_profiles_df = novel_profiles_df.append(novel_profile) #append novel profile to temporary novel profile db
				novel_profiles_df = novel_profiles_df[~novel_profiles_df.index.duplicated(keep='first')] #the current allele profile will be added to the dataframe every chunk (currently don't now how to add it after the last chunk), so need to stop build up of duplicates per cunk

		if isolate in list(novel_profiles_df.index): #if isolate have novel profile
			if len(current_st) == 0:
				current_st = profiles_chunk.iloc[-1].name #get current latest st
				current_st = re.sub('_','', str(current_st)) #remove novel st identifier (if present)
				novel_st = str(int(current_st) + 1) #next novel st
				novel_profiles_df.loc[[isolate],'ST'] = '_'+novel_st #give current novel isolate profile novel ST
				novel_sts.append(novel_st) #appened novel ST to list
				novel_isolates.append(isolate) #appened isolate with novel allele profile to list
			else:
				current_st = int(current_st) + 1 #get current latest st
				current_st = re.sub('_','', str(current_st)) #remove novel st identifier (if present) (probably needless operation, as '_' is removed first time)
				novel_st = str(int(current_st) + 1) #next novel st
				novel_profiles_df.loc[[isolate],'ST'] = '_'+novel_st #give current novel isolate profile novel ST
				novel_sts.append(novel_st) #appened novel ST to list
				novel_isolates.append(isolate) #appened isolate with novel allele profile to list
		else:
			pass

	##Variables for stats
	tot_isos = len(isolate_id)
	tot_nov = len(novel_sts)
	tot_iden = len(sts_results)
	tot_re_st = len(novel_st_multi)

	##Make table of isolates and their STs (novel or not(further down))
	run_novel_sts_out = pd.DataFrame({'Isolate_ID':novel_isolates,'ST':novel_sts})
	run_novel_sts_out.to_csv(args['out_dir']+'run_novel_sts_out.csv',index=False)

	if len(novel_st_multi) == 0:
		pass
	else:
		run_novel_multi_sts_out = pd.DataFrame({'Isolate_ID':novel_iso_multi,'ST':novel_st_multi})
		run_novel_multi_sts_out.to_csv(args['out_dir']+'run_novel_multi_sts_out.csv',index=False)

	if len(sts_results) == 0 and len(novel_st_multi):
		print('All isolates had novel STs')
	else:
		run_ident_sts_out = pd.DataFrame({'Isolate_ID':isolate_results,'ST':sts_results})
		run_ident_sts_out.to_csv(args['out_dir']+'run_ident_sts_out.csv',index=False)
		print(f'Of {tot_isos} isolate(s), {tot_nov} had novel ST(s), {tot_iden} had existing ST(s) and {tot_re_st} isolate(s) were identical in this run')

	##Save new reference db
	timestr = time.strftime("%d_%m_%Y-%H_%M")
	novel_profiles_df.set_index('ST', inplace=True)

	for index, row in novel_profiles_df.iterrows():
		a = row.to_list()
		b = index
		a.insert(0,b)
		with open(args['db'], 'a') as f:
			writer = csv.writer(f)
			writer.writerow(a)
			f.close()

	##copy file to new database name
	shutil.copy2(args['db'], f"{args['db_dir']}profiles.list.{timestr}.csv")


	# End timer
	end_time = time.time()

	# Calculate elapsed time
	elapsed_time = end_time - start_time
	print("Elapsed time: ", elapsed_time)

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

	##Create memory exceptions
	memory_limit_half()
	try:
		main(vars(args))
		sys.exit(0)
	except MemoryError:
		sys.stderr.write('\n\nERROR: Memory Exception\n')
		sys.exit(1)
