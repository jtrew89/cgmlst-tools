# This rule gets ST for isolates in run, either attributing an existing ST or making a novel on. It updates the ST database with novel STs and outputs STs of isolates in run
# Libraries used
import os
import time

# Variables/directories used
scripts_dir = "scripts"
db_dir = "dbs/"
out_dir = "etoki_out/"
timestr = time.strftime("%d_%m_%Y-%H_%M")

# Run allele profile comparison to get cgMLST
rule get_sts:
	input:
		alle_pros = os.path.join(out_dir, "results_alleles.tsv"),
		st_db = os.path.join(db_dir, "profiles.list.15_05_2026-15_25.csv")
		#st_db = rules.get_previous.output.out_sts # Supposedly uses the last db output by the rule, so something to use for later when pipeline is properly being used
	output:
		out_sts = os.path.join(db_dir,"profiles.list." + timestr + ".csv")
	
	log:
		os.path.join(out_dir, "log/get_sts" + timestr + ".log")

	shell:
		"{scripts_dir}/cg_mlst_get_st_mem.py "
		"-i {input.alle_pros} "
		"-db {input.st_db} "
		"-dd {output.out_sts} "
		"-od {out_dir} "
		"> {log} 2>&1"