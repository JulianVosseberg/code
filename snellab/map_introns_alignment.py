#!/usr/bin/python3

# Author: Sjoerd Gremmen
# Adapted by Julian Vosseberg

import sys; import getopt

def usage():
    print("\tUsage: map_introns_alignment.py [-h] [ options ] <cds_euk_ids.tsv> <alignment.fa>")
    print("\t  cds_euk_ids.tsv is the output provided by get_cds_gff_euk_id.sh\n")
    print("Options:\n\t-g: supergroups.tsv (ABBR\\tsupergroup (default: ~/julian2/snel-clan-genomes/eukarya/5_supergroups.tsv))")
    print("\t-o: output directory (default: current directory)")
    print("\t-s: intron shift allowed (default: 0)")

# Parse arguments
optlist, args = getopt.getopt(sys.argv[1:], 'g:o:s:h')
opts = {}
for k,v in optlist:
    if k == '-h':
        usage(); sys.exit()
    else:
        opts[k] = v
if len(args) != 2:
    print("Please provide a cds file and alignment\n")
    usage(); sys.exit()
else:
    information_about_CDS = args[0]
    Your_protein_alignment = args[1]

if '-g' in opts.keys():
    Species_to_supergroups = opts['-g']
else:
    Species_to_supergroups = '/hosts/linuxhome/memory/julian2/snel-clan-genomes/eukarya/5_supergroups.tsv'
if '-o' in opts.keys():
    output_path = opts['-o'] + "/"
else:
    output_path = ''
if '-s' in opts.keys():
    shifts = int(opts['-s'])
    if shifts > 10:
        print("Warning: the number of shifts you want to tolerate seems very high", file = sys.stderr)
else:
    shifts = 0

# Set thresholds for considering an intron a LECA intron
threshold_percentage_genes = 7.5 # Intron has to be present in at least this many genes
threshold_percentage_species = 15 # Intron has to be present in at least this many species
threshold_species = 2 # Intron has to be present in at least X Opimoda (~unikont) and X Diphoda (~bikont) species

# Step 1: Parse input files
# Parse CDS information; euk_cds_dict contains coordinates in format {seq_id:[['+/-', start, end], ['+/-', start, end]]}
cds_information = open(information_about_CDS,"r")    
euk_cds_dict = {}
for line in cds_information:
    line = line.strip()
    words = line.split("\t")
    sequence_id = words[0]
    start_cds = int(words[1])
    stop_cds = int(words[2])
    direction = words[3]
    if sequence_id in euk_cds_dict:
        temporary_list = [direction, start_cds,stop_cds]
        euk_cds_dict[sequence_id].append(temporary_list)
    else:
        euk_cds_dict[sequence_id]= [[direction,start_cds,stop_cds]]
cds_information.close()

# Check if all coding parts of one gene lie in the same direction
for key in euk_cds_dict:
    direction = 'none'
    for coding_sequence in euk_cds_dict[key]:
        if direction == 'none':
            direction = coding_sequence[0]
        else:
            if coding_sequence[0] != direction:
                print ("FATAL ERROR: the CDSs of ",key," are not in the same direction", file = sys.stderr)
                sys.exit()

# Parse supergroup information
supergroups_doc = open(Species_to_supergroups, "r")
supergroups ={}
groups = []
for line in supergroups_doc:
    line = line.strip()
    words = line.split()
    supergroups[words[0]] =str(words[1])
    if words[1] not in groups:
        groups.append(words[1])
supergroups_doc.close()
                
# Parse alignment
align = open(Your_protein_alignment, 'r')
aligned_proteins = {} 
alignment = ""
OG = '' 
OG_list = []
OGdict = {}
number_species_dict ={}
total_number_of_genes = {}
species_list = {}
number_of_species ={}
for line in align:
    line = line.strip()
    if line.startswith(">"): 
        end = line.rfind("_")
        OG = line[1:end]
        species = line[end+1:end+5]
        if OG not in OG_list:
            OG_list.append(OG)
            total_number_of_genes[OG] = 0
            OGdict[OG] = {}
            number_species_dict[OG] = {}
            number_of_species[OG] = 0
            species_list[OG] = []
            for group in groups:
                OGdict[OG][group]=0
                number_species_dict[OG][group] = 0
        if species not in supergroups:
            print("Error: ",line[end+1:end+5]," not in supergroup file", file = sys.stderr)
        else:
            OGdict[OG][supergroups[line[end+1:end+5]]] += 1
            total_number_of_genes[OG] +=1
            if species not in species_list[OG]:
                species_list[OG].append(species)
                number_of_species[OG] += 1
                number_species_dict[OG][supergroups[line[end+1:end+5]]]+=1
        if alignment != '':
            aligned_proteins[protein] = alignment
            alignment = ''
        protein = line[1:]  
    else:
        alignment += line
aligned_proteins[protein] = alignment 
align.close()

# Step 2: Map intron positions onto the proteins
id_length_introns = {}
for gene in euk_cds_dict: 
    id_length_introns[gene] = []
    length_without_introns = 0
    location_introns = []
    startcds = euk_cds_dict[gene][0][1]
    stopcds = euk_cds_dict[gene][len(euk_cds_dict[gene])-1][2] #-1 needed for counting in python
    # Some of the minus directed CDSs needed to be reversed.
    if euk_cds_dict[gene][0][0] == "-":
        if startcds < stopcds:
            euk_cds_dict[gene] = euk_cds_dict[gene][::-1]
    for intron_information in euk_cds_dict[gene]:
        if length_without_introns == 0:
            length_without_introns = intron_information[2] - intron_information[1] +1 
            phase = length_without_introns%3 
            location_intron = (length_without_introns-phase)/3 +1 
        else:
            intron = [phase, location_intron]
            location_introns.append(intron)
            # Adding intron to list to prevent end of protein being seen as intron.
            length_without_introns += (intron_information[2]-intron_information[1])+1 
            phase = length_without_introns%3 
            location_intron = (length_without_introns-phase)/3 +1 
        id_length_introns[gene] = [length_without_introns/3, location_introns]
        # everything is devided by three, from length of mRNA (nucleotides) to polypeptide length (amino acid). 
    if length_without_introns%3 != 0:
        print ("ERROR: ",gene, " seems to be incorrectly annotated", file = sys.stderr) #A small check if the genes are correctly annotated

# Step 3 : Map intron positions onto the alignment
if shifts <= 3:
    aa = 1    
else:
    if shifts%3 == 0: 
        aa = shifts/3
    else:
        aa =(shifts-shifts%3)/3 +1 # 0-3 --> aa=1, 4-6 --> aa=2, 7-9 --> aa=3, ...

introns_shifts = {} 
dict_with_introns ={}
gene_supergroup_dict ={}
LECA_introns = {}
dominant_phases = {}
for OG in OG_list:
    gene_supergroup_dict[OG] = {}
    dict_with_introns[OG]={}
    for protein in aligned_proteins:
        position_intron = 0 #regulates that not every intron is checked every time.
        if OG in protein: 
            begin = protein.rfind("_")+1
            ID = protein[begin:]
            location = 0
            position = 0 
            if ID in id_length_introns:
                for character in aligned_proteins[protein]: 
                    position += 1
                    if character != "-":
                        location +=1
                    if position_intron < len(id_length_introns[ID][1]):
                        if id_length_introns[ID][1][position_intron][1] == location:
                            intron= id_length_introns[ID][1][position_intron]
                            position_intron+=1 
                            if position not in dict_with_introns[OG]:
                                dict_with_introns[OG][position]=[[],[],[]]
                            dict_with_introns[OG][position][intron[0]].append(protein)
                            
    introns_shifts[OG] = {} 
    dominant_phases[OG] = {}
    for position in dict_with_introns[OG]: 
        introns_shifts[OG][position]=[]
        dominant = dict_with_introns[OG][position][0]
        dominant_phase = 0
        if len(dict_with_introns[OG][position][1])>len(dominant):
            dominant = dict_with_introns[OG][position][1]
            dominant_phase = 1 
        if len(dict_with_introns[OG][position][2]) > len(dominant): 
            dominant = dict_with_introns[OG][position][2]
            dominant_phase = 2
        dominant_phases[OG][position]=dominant_phase
        for k in range(int(position-aa),int(position+aa)+1,1):  #k is just a variable looping through the list of numbers
            if k in dict_with_introns[OG]:
                introns_shifts[OG][position].extend(dict_with_introns[OG][k])
            else:
                introns_shifts[OG][position].extend([[],[],[]]) 
        introns_shifts[OG][position] = introns_shifts[OG][position][int(3*aa+dominant_phase-shifts):int(3*aa+dominant_phase+shifts+1)]
        for j in range(0,len(introns_shifts[OG][position]),1):
            if len(introns_shifts[OG][position][j]) > len(dominant): 
                #Assumption is that in the range defined there can only be 1 LECA intron and that if shift <3 
                #there can only be one dominant intron on 1 aa position.
                del introns_shifts[OG][position]
                #Selection: If near this position is a location with more introns, this location is deleted.
                break
    for dominant_position in introns_shifts[OG]: 
        gene_supergroup_dict[OG][dominant_position] = {"Amoebozoa" : [],"Obazoa" : [], "Excavata": [], "Archaeplastida": [], "RASH": []}
        frames = -1
        for phase in introns_shifts[OG][dominant_position]:
            for supergroep in gene_supergroup_dict[OG][dominant_position]:
                gene_supergroup_dict[OG][dominant_position][supergroep].append([])
            frames += 1
            for ID in phase:
                begin = ID.rfind("_")
                species = ID[begin+1:begin+5]
                if species in supergroups:
                    try:
                        gene_supergroup_dict[OG][dominant_position][supergroups[species]][frames].append(ID)
                    except:
                        print("ERROR: something went wrong while composing the gene supergroup dictionary",ID, file = sys.stderr)

    # Write mapped introns to output file   
    OG_file = open(output_path+"analysis_"+OG+".txt", "w")
    raw = open(output_path+"raw_data_"+OG+".txt", "w")
    
    #printing a header for all documents:
    print("Intron position/Supergroup",end="\t",file=OG_file)
    print("Intron position/Supergroup",end="\t",file=raw)
    for i in range(-shifts,shifts+1,1):
        if i == 0:
            print("dominant phase", end = "\t",file = OG_file)
            print("dominant phase", end = "\t",file = raw)
        else:
            print(i,end="\t",file=OG_file)
            print(i,end="\t",file=raw)
    print("Percentage of genes","\t", "Percentage of species", file = OG_file)
    print("Percentage of genes","\t","Percentage of species", file = raw)
    
    for intron_position in sorted(gene_supergroup_dict[OG]):
        number_of_species_intron = 0
        list_with_species_intron =[]
        genes_per_position = 0
        percentage_genes_per_group={}
        percentage_species_per_group = {}
        uniconta = 0
        biconta = 0
        genes_per_phase={}
        represented_groups = 0
        for group in gene_supergroup_dict[OG][intron_position]:
            phases = 0
            percentage_genes_per_group[group] = 0
            percentage_species_per_group[group] = 0
            genes_per_group = 0
            species_per_group = 0
            for phase in gene_supergroup_dict[OG][intron_position][group]: 
                phases+=1
                if phases not in genes_per_phase:
                    genes_per_phase[phases] = 0
                
                #Counting the Opimoda and Diphoda for required threshold.
                # Counting several variables needed to calulate percentages.
                for gene in phase:
                    if group == "Amoebozoa" or group == "Obazoa":
                        uniconta += 1
                    else:
                        biconta += 1 
                    begin = gene.rfind("_")+1
                    species = gene[begin:begin+4]
                    if species not in list_with_species_intron:
                        list_with_species_intron.append(species)
                        species_per_group +=1
                        number_of_species_intron +=1
                    genes_per_phase[phases] +=1
                    genes_per_group += 1
                    genes_per_position += 1
            
            #Caluclating percentages to see if introns reach the threshold
            if genes_per_group != 0:
                percentage_genes_per_group[group] = str((genes_per_group/OGdict[OG][group])*100)[0:4]
                percentage_species_per_group[group] = str((species_per_group/number_species_dict[OG][group]*100))[0:4]
        total_percentage_of_gene = str((genes_per_position/total_number_of_genes[OG])*100)[0:4]
        percentage_of_species = str((number_of_species_intron/number_of_species[OG])*100)[0:4]
        
        #Check if introns reach the thresholds to be considered LECA genes
        if uniconta >= threshold_species and biconta >= threshold_species:
            if float(total_percentage_of_gene) > threshold_percentage_genes:
                if float(percentage_of_species) > threshold_percentage_species:
                    
                    #making the analysis file
                    k = 0
                    print(">",intron_position,"\t","Dominant_phase = ", dominant_phases[OG][intron_position], file = OG_file)
                    if intron_position in LECA_introns:
                        LECA_introns[intron_position].append(OG)
                    else:
                        LECA_introns[intron_position]=[OG]
                    for group in gene_supergroup_dict[OG][intron_position]:
                        print(group, end="\t", file = OG_file)
                        for phase in gene_supergroup_dict[OG][intron_position][group]:
                            print(phase,end ="\t",file = OG_file)
                        print(percentage_genes_per_group[group],"\t",percentage_species_per_group[group],"%",file = OG_file)
                    for phase in genes_per_phase:
                        k += 1
                        percentage_per_phase = str(genes_per_phase[phase]/genes_per_position * 100)[0:4]
                        if k == 1:
                            print("total percentage of genes", "\t",percentage_per_phase,"%",end="\t", file = OG_file)
                        else:
                            print(percentage_per_phase,"%",end="\t",file = OG_file)
                    print(total_percentage_of_gene,"%","\t", percentage_of_species,"%","\n", file = OG_file)
        print(">",intron_position,"\t","Dominant phase = ",dominant_phases[OG][intron_position], file = raw)
        
        #Builing a file with raw data
        for group in gene_supergroup_dict[OG][intron_position]:
            print(group, end="\t", file =raw)
            for phase in gene_supergroup_dict[OG][intron_position][group]:
                print(phase,end ="\t",file=raw)
            print(percentage_genes_per_group[group],file = raw)
        k = 0
        for phase in genes_per_phase:
            if genes_per_position != 0:
                k += 1
                percentage_per_phase = str(genes_per_phase[phase]/genes_per_position * 100)[0:4]
                if k > 1:
                    print(percentage_per_phase,"%",end="\t",file = raw)
                else:
                    print("total percentage of genes", "\t",percentage_per_phase,"%",end="\t", file =raw)
        print(total_percentage_of_gene,"%","\n", file = raw)         
    raw.close()    
    OG_file.close()

# Make a log file
logfile = open(output_path+"logfile"+str(OG_list)+".txt","w")
print("threshold percentage genes = ",threshold_percentage_genes,"threshold percentage species: ",threshold_percentage_species,"\t","Shifts = ", shifts, file=logfile)
for intron_position in LECA_introns:
    print("\nWe found a LECA intron on position ",intron_position," of the alignment", file= logfile)
    for OG in LECA_introns[intron_position]:
        print("This LECA intron is found in ", OG,"in Phase ",dominant_phases[OG][intron_position], file = logfile)
print("\nScore:", file = logfile)
number_of_shared_introns = {}
for OG in OG_list:
    score = {}
    number_of_shared_introns[OG]={}
    for LECA_intron in LECA_introns:
        if OG in LECA_introns[LECA_intron]:
            #Because shifts should be taken in account
            for possible_shift in range(LECA_intron-aa,LECA_intron+aa+1,1):
                if possible_shift != LECA_intron:
                    if possible_shift in LECA_introns:
                        for other_OG in LECA_introns[possible_shift]:
                            k=0
                            for hypothetical_phase in range(LECA_intron*3+dominant_phases[OG][LECA_intron]-shifts,LECA_intron*3+dominant_phases[OG][LECA_intron]+shifts+1,1):
                                if k == 0:
                                    if hypothetical_phase in range(possible_shift*3+dominant_phases[other_OG][possible_shift]-shifts,possible_shift*3+dominant_phases[other_OG][possible_shift]+shifts+1,1):
                                        k+=1
                                        if other_OG in number_of_shared_introns[OG]:
                                            number_of_shared_introns[OG][other_OG]+=1
                                            print(OG,other_OG)
                                        else:
                                            number_of_shared_introns[OG][other_OG] = 1
                else:#I deliberately used the if else statements below (the same as above), because I thought using this else would save the computer a lot of calculating
                    for other_OG in LECA_introns[LECA_intron]:
                        if other_OG in number_of_shared_introns[OG]:
                            number_of_shared_introns[OG][other_OG] += 1
                        else:
                            number_of_shared_introns[OG][other_OG] = 1
logfile.close()
table = open(output_path+"table"+str(OG_list)+".txt","w")
new_list = []
lines = {}
first_line = []
for OG in number_of_shared_introns:
    first_line.append(OG)
    lines[OG] = [OG]
    new_list.append(OG)
    for other_OG in new_list:
        if other_OG == OG:
            lines[OG].append("1.0")
        elif other_OG in number_of_shared_introns[OG]:
            #Calculating the score
            score = str(number_of_shared_introns[OG][other_OG]/(number_of_shared_introns[other_OG][other_OG]))[0:4]
            lines[OG].append(float(score))
        else:
            lines[OG].append(0.0)
k = 0
for OG in first_line:
    k+=1
    if k < len(first_line):
        print("\t", end= OG ,file=table)
    else:
        print("\t", OG, file = table)
for OG in lines:
    for value in lines[OG]:
        if value == "1.0":
            print(value,file = table)
        else:
            print(value,end ="\t", file = table)
table.close()
