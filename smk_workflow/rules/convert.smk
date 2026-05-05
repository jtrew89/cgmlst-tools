# This rule parses the etoki output in a dir to a dataframe and updates the md5SUM convert table with mdtsums from novel alleles
import os

# Variables/directories used
scripts_dir = "scripts"
sch_dir = "schemas"
out_dir = "etoki_out/"

# Get fastas
rule all:
	input:
		os.path.join(out_dir,"results_alleles.tsv")

# Generate results_alleles from 

rule etoki_convert:
	input:
		all_con = os.path.join(sch_dir,"all_convert.tab")
	output:
		parsed_etoki = os.path.join(out_dir,"results_alleles.tsv")
	shell:
		"{scripts_dir}/cg_mlst_etoki_convert.py "
		"-id {out_dir} "
		"-ct {input.all_con} "
		"-od {output.parsed_etoki}"