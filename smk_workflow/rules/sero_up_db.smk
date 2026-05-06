# Thus rule adds updates serotype specific database with latest allele profiles
# Libraries used
import os
import time

# Variables/directories used
scripts_dir = "scripts"
db_dir = "dbs/typhimur/"
out_dir = "etoki_out/"
date = time.strftime("%d_%m_%Y")

# Update serotype db with the runs latest allele profiles
rule all:
	input:
		os.path.join(db_dir,"curr_db_" + date + ".tsv")

rule sero_db_update:
	input:
		sero_db = os.path.join(db_dir, "base.tsv")
	output:
		updated_db = os.path.join(db_dir,"curr_db_" + date + ".tsv")
	params:
		out_dir = out_dir
	shell:
		"{scripts_dir}/cg_mlst_add_db.py "
		"-id {out_dir} "
		"-sdb {input.sero_db} "
		"-o {output.updated_db}"