# This rule runs EToKi.py MLSType
# Libraries used
from snakemake.io import glob_wildcards
import os

# Variables/directories used
etoki = "/home/jahcubtrew/packages/EToKi/EToKi.py"
reads_dir = "reads"
out_assem = "assemblies"

# Get fastqs
wildcards_fq = glob_wildcards(os.path.join(reads_dir,'{sample}_R1.fastq.gz'))

# Generate EToKi MLST call output
rule assembly:
	input:
		for_read = os.path.join(reads_dir,'{sample}_R1.fastq.gz'),
		rev_read = os.path.join(reads_dir,'{sample}_R2.fastq.gz')
	output:
		results_alleles = os.path.join(out_assem, "{sample}.result.fasta")
	resources:
		pair_jobs=2
	params:
		sample=lambda wc: wc.sample
	shell:
		"python {etoki} assemble "
		"--pe {input.for_read},{input.rev_read} "
		"-p {out_assem}/{params.sample}"