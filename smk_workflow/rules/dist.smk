# Rule to create distance matrix for input isolates 
# Libraries used
import os
import pandas as pd
#import time

# Variables/directories used
out_dir = "etoki_out/"
pairwise_out = "etoki_out/pairwise/"
scripts_dir = "scripts"
#timestr = time.strftime("%d_%m_%Y-%H_%M")

# Get isolate list to run distance script per isolate
result_df = pd.read_table(out_dir + "results_alleles.tsv")
isolate_ids = result_df["Isolate_ID"].unique().tolist()

# File snakemake will look for, and run if not present
rule pair_dist:
	input:
		alle_pros = os.path.join(out_dir,"results_alleles.tsv")
	output:
		dist_out = os.path.join(pairwise_out,"{isolate}.tsv")
	shell:
		"{scripts_dir}/cg_mlst_dist.py "
		"-pdb {input.alle_pros} "
		"-id {wildcards.isolate} "
		"-od {pairwise_out}/ "
		"-m"

rule combine_matrix:
	input:
		expand(os.path.join(pairwise_out, "{isolate}.tsv"), isolate=isolate_ids)
	output:
		#dist_out = os.path.join(out_dir,"allele_dist_miss" + '_' + timestr + ".tsv")
		dist_out = os.path.join(out_dir,"allele_dist_miss.tsv")
	run:
		import pandas as pd
		df = pd.concat([pd.read_csv(f, sep="\t") for f in input])
		df.to_csv(output.dist_out, sep="\t", index=False)