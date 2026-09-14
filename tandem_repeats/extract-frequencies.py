#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract tandem repeat information from GTF
"""
import glob
import pandas as pd
# Set output folder
output_folder = "~/Documents/manchester/tr/results/overview/"
# Generate a cluster assignment key
metadata_filename = "~/Documents/manchester/global/data/metadata/microreact.shared.csv"
metadata = pd.read_csv(metadata_filename, header=None, index_col=0)

# Investigate gff files
gtf_files = glob.glob("/Users/user/Documents/manchester/global/data/trf/gff/*.trf.gff3")


homopolymer_results = []
tandem_repeat_results = []
cons_size_results = []
repeat_length_results = []
homopolymer_count = []
for gtf_filename in gtf_files:
    gtf = open(gtf_filename).readlines()
    
    
    isolate = gtf[0].split("_")[0]
    if isolate not in metadata.index:
        continue
    else:
        cluster = "".join(["Cluster ",str(metadata.loc[isolate][11])])
        tr_cons_size = {}
        tr_repeat_length = {}
        homopolymer_overview = {}
        for line in gtf:
            line_split = line.split(";")
        
            # Filters:
            # Percentage match >= 85
            # Consensus sequence mutliple bases else homopolymer
            # Copy number round down
            # Use consensus size number
            
            perc_match = int(line_split[4].split("=")[1])
            if perc_match >= 85:
                cons_seq =  line_split[8].split("=")[1]
                cons_size = int(line_split[3].split("=")[1])
                copies = int(line_split[2].split("=")[1].split(".")[0])
                repeat_length = copies * cons_size
                bases = list(set(cons_seq))
                if len(bases) == 1:
                    homopolymer_results += [[isolate, cluster, bases[0], repeat_length]]
                    homopolymer_overview[repeat_length] = homopolymer_overview.get(repeat_length, 0) + 1
                else:
                    tandem_repeat_results += [[isolate, cluster, cons_size, copies, repeat_length]]
                    tr_cons_size[cons_size] = tr_cons_size.get(cons_size, 0) + 1
                    tr_repeat_length[repeat_length] = tr_repeat_length.get(repeat_length, 0) + 1
        # Convert dictionary of frequencies to main list
        for size, count in tr_cons_size.items():
            cons_size_results += [[isolate,cluster, size, count]]
        for size, count in tr_repeat_length.items():
            repeat_length_results += [[isolate, cluster, size, count]]
        for size, count in homopolymer_overview.items():
            homopolymer_count += [[isolate, cluster, size, count]]
        
            
# Create dataframe of frequency stats
consensus_size_frequency = pd.DataFrame(cons_size_results)
consensus_size_frequency.columns = ["Isolate", "Cluster", "Consensus Size", "Count"]
consensus_size_frequency.to_csv("".join([output_folder, "consensus_size_freq.csv"]), index=False)


repeat_length_frequency = pd.DataFrame(repeat_length_results)
repeat_length_frequency.columns = ["Isolate", "Cluster", "Total Length", "Count"]
repeat_length_frequency.to_csv("".join([output_folder, "repeat_length_freq.csv"]), index=False)


homopolymer_frequency = pd.DataFrame(homopolymer_count)
homopolymer_frequency.columns = ["Isolate", "Cluster", "Homopolymer Run Length", "Count"]
homopolymer_frequency.to_csv("".join([output_folder, "homopolymer_freq.csv"]), index=False)


# Create dataframe of total data
total_tandem_repeat = pd.DataFrame(tandem_repeat_results)
total_tandem_repeat.columns = ["Isolate", "Cluster", "Consensus Size", "Copies", "Repeat Length"]
total_tandem_repeat.to_csv("".join([output_folder, "total_tandem_repeat.csv"]), index=False)


total_homopolymer = pd.DataFrame(homopolymer_results)
total_homopolymer.columns = ["Isolate", "Cluster", "Nucleotide", "Repeat Length"]
total_homopolymer.to_csv("".join([output_folder, "total_homopolymer.csv"]), index=False)








