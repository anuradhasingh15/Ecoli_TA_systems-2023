#This script has been adapted from https://github.com/ghoresh11/kpneumoniae_TAs

import os


## for the prediction of toxin sequences:-
## given all the itoxin sequences
## get all of them -> create a blastn DB that is all the toxins
## using the same BLAST cutoffs used in SLING
## blast the genomes against the toxins DB and see if I find any


## step 1 -> aggregate all antitoxin sequences
toxins_dir = "/home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/sling/output/ID_1_GROUP/hits_clusters/"
toxins_files = os.listdir(toxins_dir)
#print(toxins_files)

toxins = open("all_toxins_set1.fasta", "w")
for f in toxins_files:
    if not f.endswith(".txt"):
        continue
    name = f.replace(".txt","")
    #print(name)
    line_num = 1
    with open(os.path.join(toxins_dir,f)) as f_open:
        for line in f_open:
            if line.startswith("#") or line.startswith("Strain"):
                continue
            toks = line.strip().split(",")
            #print(toks)
            toxins.write(">" + name + "_" + str(line_num) + "\n" + toks[-3] + "\n")
            line_num += 1
            
toxins.close()
