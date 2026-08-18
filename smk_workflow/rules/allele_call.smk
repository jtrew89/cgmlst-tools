# This rule runs EToKi.py MLSType
# Libraries used
from snakemake.io import glob_wildcards
import os

# Variables/directories used
etoki_dup = "/home/jahcubtrew/packages/etoki_duplicate/EToKi.py"
out_assem = "assemblies"
sch_dir = "schemas"
out_dir = "etoki_out"

# Get fastas
wildcards_fa = glob_wildcards(os.path.join(out_assem,'{sample}.result.fasta'))

# Generate EToKi MLST call output
rule allele_call:
	input:
		iso_assem = os.path.join(out_assem,'{sample}.result.fasta'),
		ref = os.path.join(sch_dir,'all_references.fasta'),
		all_con = os.path.join(sch_dir,'all_convert.tab')
	output:
		results_alleles = os.path.join(out_dir, "{sample}_results_alleles.fasta")
	resources:
		pair_jobs=3
	shell:
		"python {etoki_dup} MLSType "
		"-i {input.iso_assem} "
		"-r {input.ref} "
		"-k G749 "
		"-o {output.results_alleles} "
		"-d {input.all_con} -l 44"

# Generate file to indicate allele_call is complete
rule allele_call_complete:
	input:
		expand(os.path.join(out_dir, "{sample}_results_alleles.fasta"),sample=wildcards_fa.sample)
	output:
		complete = touch("etoki_out/.allele_call_complete")