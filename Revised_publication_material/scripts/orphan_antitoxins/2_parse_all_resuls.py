#This script has been adapted from https://github.com/ghoresh11/kpneumoniae_TAs

import os


# Set the directory paths for the input files
group_dir = "/home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/sling/output/ID_1_GROUP/partners_clusters/"
blast_dir = "/home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/orphan_antitoxins_1/blast_output/"

#my blast output format is "6 qseqid sseqid pident length qstart qend sstart send sstrand evalue bitscore'

# Define two dictionaries to store the antitoxin and strain data
antitoxins = {}
strains = {}

# Get a list of all the antitoxin cluster files
antitoxin_clusters = os.listdir(group_dir)
#print(antitoxin_clusters)

# Iterate over each antitoxin cluster file to extract data and store it in the dictionaries
for file in antitoxin_clusters:
    # Check if the file is a text file, if not skip to the next file
    if not file.endswith(".txt"):
        continue
    # Get the cluster name by removing the ".txt" extension from the file name
    cluster = file.replace(".txt", "")
    #print(cluster)
    # Add a new entry to the antitoxins dictionary for this cluster
    antitoxins[cluster] = {"hits": set(), "avg_length": 0, "domains": [], "strains": [], "new_strains": []}
    #print(antitoxins[cluster])
    
    # Open the antitoxin cluster file
    with open(os.path.join(group_dir, file)) as f_open:
        # Iterate over each line in the file
        for line in f_open:
            # Check if the line contains the domains information, if so store it in the antitoxin dictionary
            if "Domains:" in line:
                line = line.strip().split(":")[-1]
                #print(line)
                antitoxins[cluster]["domains"] = line.strip().split(",")
                #print(antitoxins[cluster]["domains"])
                continue
                # Check if the line contains the average length information, if so store it in the antitoxin dictionary
            if "Partner Length" in line:
                line = line.strip().split(":")[-1]
                #print(line)
                antitoxins[cluster]["avg_length"] = float(line)
                #print(antitoxins[cluster])
                continue
            # Check if the line is a header line or a comment line, if so skip to the next line
            if line.startswith("Strain") or line.startswith("##") or line.startswith("#"):
                continue
            # Otherwise, extract the strain and hit information from the line and store it in the antitoxin dictionary
            toks = line.strip().split(",")
            #print(toks)
            antitoxins[cluster]["hits"].add(toks[1])
            antitoxins[cluster]["strains"].append(toks[0])
            #print(antitoxins[cluster])
            
            # Add a new entry to the strains dictionary for this strain, if it does not exist already
            if toks[0] not in strains:
                strains[toks[0]] = {"new_antitoxins": []}
            # Add a new entry to the cluster dictionary for this strain, if it does not exist already
            if cluster not in strains[toks[0]]:
                strains[toks[0]][cluster] = {"contigs": [], "lengths": [], "starts": [], "stops": [], "strand": []}
                
            # Store the contig, start, stop, and strand information for this hit in the strain's cluster dictionary
            strains[toks[0]][cluster]["contigs"].append(toks[5])
            strains[toks[0]][cluster]["starts"].append(int(float(toks[9])))
            strains[toks[0]][cluster]["stops"].append(int(float(toks[10])))
            strains[toks[0]][cluster]["strand"].append(toks[6])
            
#Get a list of all the blast result files
blast_results = os.listdir(blast_dir)
print(blast_results)


def add_cluster(antitoxins, cluster, strain, s, contig, start, stop, strand):

    antitoxins[cluster]["new_strains"].append(s)
    strain["new_antitoxins"].append(cluster+"(" + contig + ":" +strand + ":" + str(start) + ":" + str(stop) + ")")
    if cluster not in strain:
        strain[cluster] = {"contigs": [], "lengths": [], "starts":[], "stops":[], "strand":[]}
    
    strain[cluster]["contigs"].append(contig)
    strain[cluster]["starts"].append(start)
    strain[cluster]["stops"].append(stop)
    strain[cluster]["strand"].append(strand)
    return
    
for file in blast_results:
    if not file.endswith(".tab"):
        continue
    #print(strains)
    
    s = file.split(".")[0]
    #print(s)
    if s in strains:
        strain = strains[s]
    # rest of the code
    else:
        print(f"Key {s} not found in strains dictionary.")

        #Key ST000120035 not found in strains dictionary because this strain doesn't contain any TA systems

    #print(strain)

    #print(strain)
    
    with open(os.path.join(blast_dir,file)) as f_open:
        for line in f_open:
            #print(line)
            #my blast output format is "6 qseqid sseqid pident length qstart qend sstart send sstrand evalue bitscore'
            toks = line.strip().split()
            #print(toks)
            cluster = toks[1].split("_")[0]
            identity = float(toks[2])
            alingment_length = int(toks[3])
            
            ## calculate the length of a protein:
            toks2 = toks[0]
            #print(toks2)
            strand = toks[8].replace("Strand:","")
            start = int(toks[4].replace("Start:",""))
            stop = int(toks[5].replace("Stop:",""))
            contig = toks[0]
            
            length = (stop - start)/3
            #print(length)
            if length < 50 or length > 150 or identity < 75:# keeping with SLING identity threshold, allow for bigger antitoxins
                continue
            
            if alingment_length <50:
                continue
            
            if cluster not in strain:
                ##found a new cluster
                add_cluster(antitoxins, cluster, strain, s, contig, start, stop, strand)
                
            else:
                if contig not in strain[cluster]["contigs"]:
                    add_cluster(antitoxins, cluster, strain, s, contig, start, stop, strand)
                else: ## same contig, check start stop and strand
                    new = True
                    indices = [i for i, x in enumerate(strain[cluster]["contigs"]) if x == contig]
                    for i in indices:
                        curr_start = strain[cluster]["starts"][i]
                        curr_stop = strain[cluster]["stops"][i]
                        curr_strand = strain[cluster]["strand"][i]
                        if curr_strand == strand and (curr_start <= start + 100 and curr_start >= start - 100) \
                        or (curr_stop <= stop + 100 and curr_stop >= stop - 100):
                            #print("Stops: %d vs %d" % (curr_stop, stop))
                            #print("Starts: %d vs %d" % (curr_start, start))
                            #print("Strands: %s vs %s" % (curr_strand, strand))
                            new = False
                            continue ## it's the same ORF
                    if new:
                        add_cluster(antitoxins, cluster, strain, s, contig, start, stop, strand)
                        
                        
#sftp://nasls@10.2.0.53/home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/orphan_antitoxins_1/orphan_at_results
out = open("/home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/orphan_antitoxins_1/orphan_at_results/summary_per_antitoxin.csv","w")
out.write("ID, Domains, Num_copies, Average_Length, Num_new_hits\n")
for c in antitoxins:
    out.write(",".join(map(str,[c, ";".join(antitoxins[c]["domains"]), len(antitoxins[c]["strains"]),\
    antitoxins[c]["avg_length"], len(antitoxins[c]["new_strains"]) ])) + "\n")
out.close()


out = open("/home/nasls/ANU/phd_project/TA-SYSTEMS/TA_systems_clean/set1/orphan_antitoxins_1/orphan_at_results/summary_per_strain.csv","w")
out.write("Name, Num_new_hits, IDs\n" )
for s in strains:
    out.write(",".join(map(str, [s, len(strains[s]["new_antitoxins"]), ";".join(strains[s]["new_antitoxins"])])) + "\n")

out.close()
            
            
            
            
            

