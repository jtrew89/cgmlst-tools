# This rule runs EToKi.py MLSType
from snakemake.io import glob_wildcards
import os

# Variables/directories used
etoki = "/home/jahcubtrew/packages/etoki_duplicate/EToKi.py"
genome_dir = "assemblies"
sch_dir = "schemas"
out_dir = "etoki_out"

# Get fastas
wildcards_fa = glob_wildcards(os.path.join(genome_dir,'{sample}.fasta'))

# Generate EToKi MLST call output
rule all:
	input:
		expand(os.path.join(out_dir,'{sample}_results_alleles.fasta'), sample = wildcards_fa.sample)

rule allele_call:
	input:
		iso_assem = os.path.join(genome_dir,'{sample}.fasta'),
		ref = os.path.join(sch_dir,'all_references.fasta'),
		all_con = os.path.join(sch_dir,'all_convert.tab')
	output:
		results_alleles = os.path.join(out_dir, "{sample}_results_alleles.fasta")
	shell:
		"python {etoki} MLSType "
		"-i {input.iso_assem} "
		"-r {input.ref} "
		"-k G749 "
		"-o {output.results_alleles} "
		"-d {input.all_con} -l 44"