#!/usr/bin/env python

"""
Script written to run chebbaca's allele call command together with cg_mlst_INF_list.py and
cg_mlst_insert_inf.py to list all novel alleles and mark them before inserting them into
the loci schema. This is so that novel allele can be identified in any updated schemas.
"""

##Import modules
import argparse
import cg_mlst_INF_list as inf_lst
import cg_mlst_insert_inf as ins_inf
import cg_mlst_get_st as get_st
import cg_mlst_add_db as up_db
import subprocess
import os
import time

##Create arguments
parser = argparse.ArgumentParser(description='Script written to run chebbaca allele call command together with cg_mlst_INF_list.py and cg_mlst_insert_inf.py to list all novel alleles and mark them before inserting them into the loci schema. This is so that novel allele can be identified in any updated schemas.')
parser.add_argument('-id', '--input_dir', dest='in_dir', help='Path to directory with input assemblies', required=True)
parser.add_argument('-sd', '--schema_dir', dest='sch_dir', help='Path to directory with schema loci fasta files', required=True)
parser.add_argument('-t', '--train', dest='train', help='Training file to be used with prodigal, with path', required=False)
parser.add_argument('-od', '--output_directory', dest='out_dir', help='Path to output directory', required=True)
parser.add_argument('-c', '--cpus', dest='cpus', help='CPUs to use in allele call', required=False, type=int)
parser.add_argument('-dd', '--db_dir', dest='db_dir', help='Directory to output update db', required=True)
parser.add_argument('-db', '--profiles_db', dest='db', help='Path to profiles.list (whichever version) and file name', required=True)
parser.add_argument('-sdb','--strain_db',dest='str_db',help="Strain db to update and path to",required=True)

args = parser.parse_args()

##Variables used in script
date = time.strftime("%d%m%Y")
run_id = os.path.basename(os.path.dirname(args.out_dir)) + '/'

##Run chewBBACA
with open(f'{args.dir}+{date}+.log', 'w') as log_file:
	subprocess.run(
		['chewBBACA.py', 'AlleleCall',
		'-i', args.in_dir, f'-g {args.sch_dir}',
		'--no-inferred --output-novel',
		f'--ptf {args.train}', f'-o {args.out_dir}.{run_id}',
		f'--cpu {cpus}'
		],
		stdout=log_file
			)

##Update serovar database with new allele profiles
update_db_dic = {
	'in_dir':args.out_dir,
	'str_db':args.str_db
		}
up_db.main(update_db_dic)

##Get list of novel alleles
novel_alle_dic = {
	'in_dir':args.out_dir,
	'out_dir':args.out_dir
		}
inf_lst.main(novel_alle_dic)

##Update alleles schema fasta files with new alleles marked with inf
ins_alle_dic = {
	'in_dir':args.out_dir,
	'schm_dir':args.sch_dir
		}
ins_inf.main(ins_alle_dic)

##Get STs for isolates just run and/or update profile.list database with novel ST profiles
get_st_dic = {
	'in_query':args.out_dir+'results_allales.tsv',
	'db':args.db, 'db_dir':args.db_dir,
	'out_dir':args.out_dir
		}
get_st.main(get_st_dic)
