#!/usr/bin/env python

# Authors: Sjoerd Gremmen, Michelle Schinkel and Julian Vosseberg

# Load modules
import sys
import os
import argparse
import gzip
from eukarya import supergroups5 as supergroups
import csv

# Functions
def parse_alignment(fasta_file):
    species_seqids_dict = {}
    OG_dict = {}
    aligned_proteins = {}
    alignment = ""
    for line in fasta_file:
        line = line.rstrip()
        if line[0] == '>':
            if alignment != '':
                aligned_proteins[seqid] = alignment
                alignment = ''
            parts = line.split('_')
            OG = parts[0][1:]
            seqid = parts[1]
            species = seqid[:4]
            species_seqids_dict[species] = species_seqids_dict.get(species, []) + [seqid]
            OG_dict[OG] = OG_dict.get(OG, []) + [seqid]
            if species not in supergroups:
                sys.stderr.write(f"Warning: {species} not recognised.\n")
        else:
            alignment += line
    aligned_proteins[seqid] = alignment
    return species_seqids_dict, OG_dict, aligned_proteins

def get_coordinates(species_seqids_dict, euk_path = '/home/julian/julian2/snel-clan-genomes/eukarya_new'):
    seqid_coordinates = {}
    for species, seqids in species_seqids_dict.items():
        if species in ('BBRI', 'ESIL'):
            sys.stderr.write(f"Warning: no CDS coordinates can be retrieved from the gff file for {species}.\n")
            continue
        if os.path.isfile(f'{euk_path}/data_set/gff_files/{species}.gff3.gz'):
            gff_file_name = f'{euk_path}/data_set/gff_files/{species}.gff3.gz'
        elif os.path.isfile(f'{euk_path}/data_set/gff_files/{species}.gff.gz'):
            gff_file_name = f'{euk_path}/data_set/gff_files/{species}.gff.gz'
        else:
            sys.stderr.write(f'Warning: no GFF file found for {species}.\n')
            continue
        original_ids = {}
        print(f'Reading metadata file for {species}...')
        with open(f'{euk_path}/data_set/proteomes_metadata/{species}.metadata.txt') as metadata_file:
            metadata_file.readline()
            lines = csv.reader(metadata_file, delimiter = '\t')
            for fields in lines:
                if fields[0] in seqids:
                    original_id = fields[8]
                    if original_id[:3] == 'ENS':
                        original_id = original_id[:original_id.find('.')]
                    if species == 'ODIO':
                        original_id = 'GSOID_T' + original_id[6:]
                    elif species == 'BSCH':
                        original_id = original_id.replace('_', '.t')
                    elif species == 'PMIN':
                        original_id = original_id.replace('-RA', '-tr')
                    elif species == 'AGAM':
                        original_id = original_id.replace('-RA', '-PA')
                    elif species == 'KFLA':
                        original_id += '.CDS'
                    elif species in ('LCOR', 'MVER'):
                        original_id = 'mRNA_' + original_id
                    elif species in ('CVEL', 'GINT'):
                        original_id = original_id[:-3] # Remove -p1
                    elif species in ('SFAL', 'ACOE', 'VCAR'):
                        original_id = original_id[:-2] # Remove .p
                    elif species in ('SMOE', 'CSUB', 'MSPE'):
                        original_id = fields[10]
                    original_ids[original_id] = fields[0]
                if len(original_ids) == len(set(seqids)):
                    break
            else:
                sys.stderr.write(f'Warning: not all sequence IDs for {species} found in metadata file.\n')
        if species in ('SPLU', 'CANG', 'MLAR', 'LTRA', 'ODIO'):
            with gzip.open(gff_file_name) as gff_file:
                for line in gff_file:
                    line = line.decode().rstrip()
                    if line[0] == '#':
                        continue
                    fields = line.split('\t')
                    if species == 'ODIO' and fields[2] == 'gene' or species != 'ODIO' and fields[2] == 'mRNA':
                        info = fields[8].split(';')
                        for field in info:
                            key, value = field.split('=')
                            if key == 'ID':
                                new_id = value
                                if species == 'ODIO':
                                    new_id = 'rna' + new_id[4:] # gene --> rna
                            elif species == 'ODIO' and key == 'Name' or species != 'ODIO' and key == 'proteinId':
                                try:
                                    eukarya_id = original_ids.pop(value)
                                    original_ids[new_id] = eukarya_id
                                except KeyError:
                                    pass
                                break
        print(f'Reading gff file for {species}...')
        with gzip.open(gff_file_name) as gff_file:
            version3 = True
            if gff_file_name.endswith('gff.gz'):
                version3 = False
            for line in gff_file:
                line = line.decode()
                if line[0] == '#':
                    continue
                line = line.rstrip().rstrip(';')
                if species == 'BSCH' and line == '':
                    continue
                fields = line.split('\t')
                if fields[2] == 'CDS':
                    if version3: #or species in ('NGAD', 'KFLA', 'GINT'):
                        info = fields[8].split(';')
                        sep = '='
                    else:
                        info = fields[8].split('; ')
                        sep = ' '
                    info_dict = {}
                    for field in info:
                        if sep not in field:
                            #sys.stderr.write(f'Warning: something odd with separator ({sep}) in this line in the gff file for {species}:\n{line}\n')
                            continue
                        key, value = field.split(sep)
                        info_dict[key] = value
                    if 'protein_id' in info_dict and species != 'ODIO':
                        original_id = info_dict['protein_id']
                    elif 'proteinId' in info_dict:
                        original_id = info_dict['proteinId']
                    elif 'Parent' in info_dict and species != 'OVOL':
                        original_id = info_dict['Parent']
                    elif 'transcript_id' in info_dict:
                        original_id = info_dict['transcript_id'][1:-1] # Remove ""
                    elif 'ID' in info_dict:
                        original_id = info_dict['ID']
                    else:
                        sys.stderr.write(f'Error: elements in this line in the gff file for {species} not recognised:\n{line}\n')
                        break
                    if species == 'NGAD':
                        original_id = original_id.split(' ')[0]
                    elif species in ('SMOE', 'CSUB', 'MSPE'):
                        original_id = original_id[:original_id.find('.')]
                    elif species in ('SFAL', 'ACOE', 'VCAR'):
                        original_id = original_id[:original_id.find('.v')]
                    elif species == 'NGAD':
                        original_id = original_id[:original_id.find(' ')]
                    elif species in ('OVOL', 'PPAC', 'SRAT', 'GSAL', 'SBAT', 'SMAN', 'EMUL', 'TASI', 'SMED'):
                        original_id = original_id[original_id.find(':') + 1:]
                    if original_id in original_ids:
                        eukarya_id = original_ids[original_id]
                        new_entry = [int(fields[3]), int(fields[4]), fields[6]]
                        try:
                            seqid_coordinates[eukarya_id].append(new_entry)
                        except KeyError:
                            seqid_coordinates[eukarya_id] = [new_entry]
        for seqid in seqids:
            if seqid not in seqid_coordinates:
                sys.stderr.write(f'Warning: {seqid} not detected in gff file.\n')
        print('Done!')
    return seqid_coordinates

def map_introns_aa(euk_cds_dict):
    location_introns = {}
    for seqid, CDSs in euk_cds_dict.items():
        # Check if all coding parts of one gene lie in the same direction
        directions = set([cds[-1] for cds in CDSs])
        if len(directions) != 1:
            sys.exit(f"FATAL ERROR: the CDSs of {seqid} are not in the same direction.")
        location_introns[seqid] = []
        length_without_introns = 0
        #location_introns = []
        startcds = CDSs[0][0]
        stopcds = CDSs[-1][1] #-1 needed for counting in python
        # Some of the minus directed CDSs needed to be reversed.
        if directions == {"-"}:
            if startcds < stopcds:
                #euk_cds_dict[seqid] = euk_cds_dict[seqid][::-1]
                CDSs.reverse()
        # Loop through CDSs, excluding the last to prevent the end of protein being seen as intron, and add to list.
        for intron_information in CDSs[:-1]:
            length_without_introns += intron_information[1] - intron_information[0] + 1
            phase = length_without_introns % 3
            location_intron = int((length_without_introns - phase)/3 + 1)
            location_introns[seqid].append([phase, location_intron])
            # everything is devided by three, from length of mRNA (nucleotides) to polypeptide length (amino acid).
        if length_without_introns%3 != 0:
            sys.stderr.write(f"Warning: {seqid} seems to be incorrectly annotated.\n") #A small check if the genes are correctly annotated
    return location_introns

def count_groups(seqids, unique = True):
    group_counts = {g : 0 for g in set(supergroups.values())}
    species = [seqid[:4] for seqid in seqids]
    if unique:
        species = set(species)
    for spec in species:
        group_counts[supergroups[spec]] += 1
    return group_counts

# Parse arguments
parser = argparse.ArgumentParser(description = "This script maps intron positions onto a protein alignment.")
parser.add_argument("alignment", help = 'protein alignment in fasta format')
parser.add_argument('-o', metavar = 'outdir', help = 'directory for output files (default: current)')
parser.add_argument('-e', metavar = 'eukarya_path', help = 'directory containing the Eukarya database (default: ~julian/julian2/snel-clan-genomes/eukarya)')
parser.add_argument('-i', help = 'infer LECA introns and introns predating duplications (default: off)', action = 'store_true')
parser.add_argument('-s', metavar = 'shifts', help = 'intron shift allowed (default: 0)', type = int, default = 0)
parser.add_argument("-p", metavar = 'gene%', help = "percentage of genes in which an intron at least has to occur to call it a LECA intron (default: 7.5)", default = 7.5)
parser.add_argument("-t", metavar = 'species%', help = "percentage of species in which an intron at least has to occur to call it a LECA intron (default: 15)", default = 15)
parser.add_argument("-n", metavar = 'species', help = "number of Opimoda and Diphoda species an intron at least has to occur to call it a LECA intron (default: 2)", default = 2)
#parser.add_argument('-d', metavar = 'domains', help = 'Pfam domains instead of full-length genes', action = 'store_true')
args = parser.parse_args()

# To add: sequence ID --> paralogue table

if args.o:
    output_path = args.o
    if not os.path.isdir(output_path):
        os.makedirs(output_path)
else:
    output_path = '.'
if args.e:
    euk_path = args.e
else:
    euk_path = '/home/julian/julian2/snel-clan-genomes/eukarya'
inference = False
if args.i:
    inference = True
shifts = args.s
if shifts > 5:
    sys.stderr(f"Warning: the number of shifts you want to tolerate ({shifts}) seems very high\n.")

# Set thresholds for considering an intron a LECA intron
threshold_percentage_genes = args.p # Intron has to be present in at least this many genes
threshold_percentage_species = args.t # Intron has to be present in at least this many species
threshold_species = args.n # Intron has to be present in at least X Opimoda (~unikont) and X Diphoda (~bikont) species

groups = set(supergroups.values())
# Parse alignment and create species_seqids_dict
# Obtain coordinates

# Step 1: Parse input files

# Parse alignment
with open(args.alignment) as fasta_file:
    species_seqids_dict, OG_dict, aligned_proteins = parse_alignment(fasta_file)
OG_group_seq_count = {OG : count_groups(seqs, unique = False) for OG, seqs in OG_dict.items()}
OG_group_species_count = {OG : count_groups(seqs, unique = True) for OG, seqs in OG_dict.items()}

# Get CDS coordinates
euk_cds_dict = get_coordinates(species_seqids_dict)

# Step 2: Map intron positions onto the proteins
location_introns = map_introns_aa(euk_cds_dict)

# Step 3: Map intron positions onto the alignment
if shifts <= 3:
    aa = 1
else:
    if shifts%3 == 0:
        aa = shifts/3
    else:
        aa =(shifts - shifts%3)/3 + 1 # 0-3 --> aa=1, 4-6 --> aa=2, 7-9 --> aa=3, ...

introns_shifts = {}
dict_with_introns = {} # To be renamed
gene_supergroup_dict ={}
LECA_introns = {}
dominant_phases = {}
for OG, seqids in OG_dict.items():
    gene_supergroup_dict[OG] = {}
    dict_with_introns[OG] = {}
    for seqid in seqids:
        if seqid in location_introns:
            position_intron = 0 #regulates that not every intron is checked every time.
            location = 0
            for position, character in enumerate(aligned_proteins[seqid]):
                if character == "-":
                    continue
                location += 1
                if position_intron < len(location_introns[seqid]):
                    if location_introns[seqid][position_intron][1] == location:
                        intron = location_introns[seqid][position_intron]
                        position_intron += 1
                        if position not in dict_with_introns[OG]:
                            dict_with_introns[OG][position]=[[],[],[]]
                        dict_with_introns[OG][position][intron[0]].append(seqid)

    introns_shifts[OG] = {}
    dominant_phases[OG] = {}
    for position in dict_with_introns[OG]:
        introns_shifts[OG][position] = []
        dominant = dict_with_introns[OG][position][0]
        dominant_phase = 0
        if len(dict_with_introns[OG][position][1]) > len(dominant):
            dominant = dict_with_introns[OG][position][1]
            dominant_phase = 1
        if len(dict_with_introns[OG][position][2]) > len(dominant):
            dominant = dict_with_introns[OG][position][2]
            dominant_phase = 2
        dominant_phases[OG][position] = dominant_phase
        for k in range(int(position-aa), int(position+aa) + 1, 1):  #k is just a variable looping through the list of numbers
            if k in dict_with_introns[OG]:
                introns_shifts[OG][position].extend(dict_with_introns[OG][k])
            else:
                introns_shifts[OG][position].extend([[],[],[]])
        introns_shifts[OG][position] = introns_shifts[OG][position][int(3*aa+dominant_phase-shifts):int(3*aa+dominant_phase+shifts+1)]
        for j in range(0, len(introns_shifts[OG][position]), 1):
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
            for supergroup in gene_supergroup_dict[OG][dominant_position]:
                gene_supergroup_dict[OG][dominant_position][supergroup].append([])
            frames += 1
            for ID in phase:
                begin = ID.rfind("_")
                species = ID[begin+1:begin+5]
                if species in supergroups:
                    try:
                        gene_supergroup_dict[OG][dominant_position][supergroups[species]][frames].append(ID)
                    except:
                        sys.stderr.write(f"Warning: something went wrong while composing the gene supergroup dictionary for {ID}.\n")

    # Write mapped introns to output file
    OG_file = open(f"{output_path}/analysis_{OG}.txt", "w")
    raw = open(f"{output_path}/raw_data_{OG}.txt", "w")

    #printing a header for all documents:
    print("Intron position/Supergroup", end = "\t",file = OG_file)
    print("Intron position/Supergroup", end = "\t",file = raw)
    for i in range(-shifts, shifts+1, 1):
        if i == 0:
            print("dominant phase", end = "\t",file = OG_file)
            print("dominant phase", end = "\t",file = raw)
        else:
            print(i, end = "\t", file = OG_file)
            print(i, end = "\t", file = raw)
    print("Percentage of genes", "\t", "Percentage of species", file = OG_file)
    print("Percentage of genes", "\t", "Percentage of species", file = raw)

    for intron_position in sorted(gene_supergroup_dict[OG]):
        #number_of_species_intron = 0
        list_with_species_intron =[]
        genes_per_position = 0
        percentage_genes_per_group={}
        percentage_species_per_group = {}
        opimoda = 0
        diphoda = 0
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
                        opimoda += 1
                    else:
                        diphoda += 1
                    #begin = gene.rfind("_")+1
                    #species = gene[begin:begin+4]
                    species = gene[:4]
                    if species not in list_with_species_intron:
                        list_with_species_intron.append(species)
                        species_per_group +=1
                        #number_of_species_intron +=1
                    genes_per_phase[phases] +=1
                    genes_per_group += 1
                    genes_per_position += 1

            #Caluclating percentages to see if introns reach the threshold
            if genes_per_group != 0:
                percentage_genes_per_group[group] = str((genes_per_group/OG_group_seq_count[OG][group])*100)[0:4]
                percentage_species_per_group[group] = str((species_per_group/OG_group_species_count[OG][group]*100))[0:4]
        total_percentage_of_gene = str((genes_per_position/len(OG_dict[OG]))*100)[0:4]
        percentage_of_species = str((len(list_with_species_intron)/sum(OG_group_species_count[OG].values()))*100)[0:4]

        #Check if introns reach the thresholds to be considered LECA genes
        if opimoda >= threshold_species and diphoda >= threshold_species:
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
                            print(','.join(phase),end ="\t",file = OG_file)
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
logfile = open(output_path+"/logfile.txt","w")
print("threshold percentage genes = ",threshold_percentage_genes,"threshold percentage species: ",threshold_percentage_species,"\t","Shifts = ", shifts, file=logfile)
for intron_position in LECA_introns:
    print("\nWe found a LECA intron on position ",intron_position," of the alignment", file= logfile)
    for OG in LECA_introns[intron_position]:
        print("This LECA intron is found in ", OG,"in Phase ",dominant_phases[OG][intron_position], file = logfile)
print("\nScore:", file = logfile)
number_of_shared_introns = {}
for OG in OG_dict:
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
table = open(output_path+"/table.txt","w")
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
