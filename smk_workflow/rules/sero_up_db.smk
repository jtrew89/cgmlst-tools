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
rule sero_db_update:
	input:
		sero_db = os.path.join(db_dir, "base.tsv"),
		alle_prof = rules.etoki_convert.output.parsed_etoki
	output:
		updated_db = os.path.join(db_dir,"curr_db_" + date + ".tsv")
	shell:
		"{scripts_dir}/cg_mlst_add_db.py "
		"-id {input.alle_prof} "
		"-sdb {input.sero_db} "
		"-o {output.updated_db}"