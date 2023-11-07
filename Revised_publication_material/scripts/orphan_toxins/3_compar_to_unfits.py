#This script has been adapted from https://github.com/ghoresh11/kpneumoniae_TAs

import os

# Step 1: generate a dictionary of strains and their unfits
# For each unfit, we need start, stop, strand, and reason
unfit_dir = "/home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/sling/output/ID_1_GROUP/unfit_clusters/"
unfit_files = os.listdir(unfit_dir)
#SSprint(unfit_files)

strains = {}

for f in unfit_files:
    if not f.endswith(".txt"):  
        continue
    with open(os.path.join(unfit_dir, f)) as f_open:
        for line in f_open:
            if "#  ID:" in line:
                line = line.strip().split(":")[-1]
                #print(line)
                toks_temp = line.strip().split("UNFIT")
                #print(toks_temp)
                continue
            if line.startswith("Strain") or line.startswith("#"):
                continue
            toks = line.strip().split(",")
            #print(toks)
            strain = toks[0]

            if strain not in strains:
                strains[strain] = {"unfits": [], "contigs": [], "starts": [], "stops": [], "strands": [], "reasons": []}

            strains[strain]["unfits"].append(toks_temp[0])
            strains[strain]["contigs"].append(toks[4])
            strains[strain]["starts"].append(int(float(toks[7])))
            strains[strain]["stops"].append(int(float(toks[8])))
            strains[strain]["strands"].append(toks[5])
            strains[strain]["reasons"].append(toks[10])
            #print(strains)
            
# Step 2: Go over all the extra hits of toxins in these strains
# and see if the toxin was found by one of the discarded toxins

toxins_out = {}
strains_out = {}

with open("/home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/orphan_toxins/orphan_at_results/summary_per_strain.csv") as f_open:
    for line in f_open:
        #print(line)
        if line.startswith("Name"):
            continue
        toks = line.strip().split(",")
        #print(toks)
        strain = toks[0]
        
        strains_out[toks[0]] = {"Num_with_unfit": 0, "Num_without_unfit": 0}
        hits = toks[2].split(";")
        #print(hits)
        for h in hits:
            toks2 = h.split("(")
            #print(toks2)
            ID = toks2[0]
            #print(ID)

            if ID not in toxins_out:
                toxins_out[ID] = {"Num_with_unfit": 0, "Num_without_unfit": 0, "hit length": 0, "Upstream length": 0, \
                "No adjacent upstream ORF": 0, "Downstream length": 0, "No adjacent downstream ORF": 0}

            toks3 = toks2[1].split(":")
            #print(toks3)
            contig = toks3[0]
            strand = toks3[1]
            start = int(toks3[2])
            stop = int(toks3[3].replace(")",""))
            flag = False
            print(strains[strain]["contigs"])
            print(len(strains[strain]["contigs"]))
            for i in range(0, len(strains[strain]["contigs"])):
                if contig == strains[strain]["contigs"][i] and strand == strains[strain]["strands"][i] and \
                strains[strain]["starts"][i] - 1000 >= start  and start <= strains[strain]["starts"][i] + 1000:
                   flag = True
                   break
            if flag:
                toxins_out[ID]["Num_with_unfit"] += 1
                strains_out[toks[0]]["Num_with_unfit"] += 1
                for r in strain["reasons"][i]:
                    if r == "":
                        continue
                    toxins_out[ID][r] += 1
            else:
                toxins_out[ID]["Num_without_unfit"] += 1
                #print(toxins_out[ID])
                strains_out[toks[0]]["Num_without_unfit"] += 1


# Step 3
with open("/home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/orphan_toxins/orphan_at_results/full_summary_per_toxin.csv", "w") as out:
    out.write("ID, Num_New, Num_with_unfit, Num_without_unfit, hit length, Upstream length, No adjacent upstream ORF, Downstream length, No adjacent downstream ORF\n")
    for ID in toxins_out:
        a = toxins_out[ID]
        out.write(",".join(map(str,[ID, a["Num_with_unfit"] + a["Num_without_unfit"], a["Num_with_unfit"], \
            a["Num_without_unfit"], a["hit length"], a["Upstream length"], \
            a["No adjacent upstream ORF"], a["Downstream length"], a["No adjacent downstream ORF"]] )) + "\n")

with open("/home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/orphan_toxins/orphan_at_results/full_summary_per_strain.csv", "w") as out:
    out.write("Name, Num_New, Num_with_unfit, Num_without_unfit\n")
    for ID in strains_out:
        s = strains_out[ID]
        out.write(",".join(map(str, [ID, s["Num_with_unfit"] + s["Num_without_unfit"],\
            s["Num_with_unfit"], s["Num_without_unfit"]])) + "\n")

                
