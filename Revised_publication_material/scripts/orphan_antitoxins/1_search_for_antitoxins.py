#This script has been adapted from https://github.com/ghoresh11/kpneumoniae_TAs

import os


## for the prediction of antitoxin sequences:-
## given all the antitoxin sequences
## get all of them -> create a blastn DB that is all the antitoxins
## using the same BLAST cutoffs used in SLING
## blast the genomes against the antitoxins DB and see if I find any


## step 1 -> aggregate all antitoxin sequences
antitoxins_dir = "/home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/sling/output/ID_1_GROUP/partners_clusters/"
antitoxin_files = os.listdir(antitoxins_dir)
#print(antitoxin_files)

antitoxins = open("all_antitoxins_set1.fasta", "w")
for f in antitoxin_files:
    if not f.endswith(".txt"):
        continue
    name = f.replace(".txt","")
    #print(name)
    line_num = 1
    with open(os.path.join(antitoxins_dir,f)) as f_open:
        for line in f_open:
            if line.startswith("#") or line.startswith("Strain"):
                continue
            toks = line.strip().split(",")
            antitoxins.write(">" + name + "_" + str(line_num) + "\n" + toks[-3] + "\n")
            line_num += 1
            
antitoxins.close()
