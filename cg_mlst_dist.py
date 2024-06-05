#!/usr/bin/python

##Import modules
import pandas as pd
import argparse

##load data
alle_prof_df = pd.read_table('C:/Users/m1011068/OneDrive - Defra/Documents/Publications/Infantis for Shaun/New_report_tree/all_results_alleles_form_new.tsv').set_index('FILE')

##Variables used in script
isolates = list(alle_prof_df.index)
out_alle_df = pd.DataFrame(index=isolates, columns=isolates)

##Replace acronyms from chewbac output and treat as missing data '0'
alle_prof_df.replace(
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
##loop throgh getting pairwise distance (allele difference) for all isolates
for isolate_1 in isolates:
    
    #new_isolates = isolates.remove(isolate_1)
    for isolate_2 in isolates:
        dist = len(alle_prof_df.loc[isolate_1].compare(alle_prof_df.loc[isolate_2]))
        
    #out_alle_df.append({'Isolates': isolate_1, isolate_2: dist[0]})
    out_alle_df.at[isolate_1, isolate_2] = dist
    
out_alle_df.set_index('Isolates')
