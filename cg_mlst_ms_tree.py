#!/usr/bin/env python
"""
mstree clusters will be based off of shared core genomic alleles. So, isolates
will only be comapred on shared cgmlst profiles. If loci are missing, they will be
omitted. Input will be a list of isolates to be an ms tree from. The script will
extract the allele profiles for those isolates.
"""

##Import modules
import pandas as pd
import argparse
import time
from tqdm import tqdm
import subprocess
import os
import sys

##Set arguments
parser = argparse.ArgumentParser(description='Script written to take a list of isolates and run reportree MSTREE2 with metadata file for isolates. Though the script can also be run to just pull allele profiles for isolates of interest, in that case do not add argument for mst')
parser.add_argument(
	'-pdb', '--profile_db', dest='pro_db',
	help='Alelle profile db with isolates you wish analysed (Include path)',
	required=True
	)
parser.add_argument(
	'-i', '--input_isolates', dest='in_iso',
	help="List of isolates to be analysed (single column, with 'Isolates' as header)",
	required=True
	)
parser.add_argument(
	'-od', '--output_dir', dest='out_dir',
	help='Directory to output results',
	required=True
	)
parser.add_argument(
	'-o', '--output_name', dest='out',
	help='Run name for output results',
	required=True
	)
parser.add_argument(
	'-mst', '--mst_medium', dest='mst',
	help='Path to mst analytical tool (currently, this only uses ReporTree)',
	required=False
	)
parser.add_argument(
	'-met', '--meta_data', dest='meta',
	help='Path to and name of metadata file to be used with mst. Must be tab-delimited',
	required= '-mst' in sys.argv
	)
parser.add_argument(
	'-c', '--meta_columns', dest='col',
	help="Column names for metadata you wish to be included in annotations of tree and some summary stats (See Reportree for description). Make sure column names are exactly as the appear in metadata file. If multiple columns, separate names with ',' no spaces. For example: -c Location,Host,AMR",
	required= '-mst' in sys.argv
	)

args = parser.parse_args()

##Load files used in analysis and set time/date
date = time.strftime("%d_%m_%Y-%H_%M")
if 'csv' in args.pro_db:
	alle_db = pd.read_csv(args.pro_db, low_memory=False)
elif 'tsv' in args.pro_db:
	alle_db = pd.read_table(args.pro_db, low_memory=False)
else:
	print(f'Input file for allele profile database ({args.pro_db}) must have a suffix of either "tsv" or "csv", and formated as such')
	exit()

isolates = list(pd.read_table(args.in_iso)['Isolates'])

##Variables used in script
iso_profs = pd.DataFrame()

##Replace acronyms from chewbac output and treat as missing data '0'
alle_db.replace(
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

##Pull profiles for each isolate from input isolates list
for isolate in tqdm(isolates):
	iso_prof = alle_db[alle_db['FILE'] == isolate]
	iso_profs = iso_profs.append(iso_prof,ignore_index=True)

iso_profs.to_csv(args.out_dir+args.out+'_profs.tsv',sep='\t',index=False)

##create minimum spanning tree using metadata and iso_profs df
if args.mst:

	os.chdir(args.mst)
	with open(f'{args.out_dir}{date}.log', 'w') as log_file:
		subprocess.run(
			[
			'./reportree.py', '-m', f'{args.out_dir}{args.meta}',
			'-a', f'{args.out_dir}{args.out}_profs.tsv', '-out', f'{args.out_dir}{args.out}',
			'--method', 'MSTreeV2', '--matrix-4-grapetree', '--columns_summary_report',
			f'{args.col}', '--metadata2report', f'{args.col}',
			'--frequency-matrix', f'{args.col}',
			'--analysis', 'grapetree'
			],
			stdout=log_file
				)
##
