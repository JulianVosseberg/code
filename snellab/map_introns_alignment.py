#!/usr/bin/env python

# Authors: Sjoerd Gremmen, Michelle Schinkel and Julian Vosseberg

# Load modules
import sys
import os
import argparse
import gzip
from eukarya import supergroups5 as supergroups
import csv
import copy
import time

# Functions
def parse_alignment(fasta_file):
    species_seqids_dict = {}
    OG_dict = {}
    aligned_proteins = {}
    alignment = ""
    aln_length = None
    for line in fasta_file:
        line = line.rstrip()
        if line[0] == '>':
            if alignment != '':
                aligned_proteins[seqid] = alignment
                if aln_length is None:
                    aln_length = len(alignment)
                elif len(alignment) != aln_length:
                    sys.exit('Error: not all aligned sequences have the same length! Analysis aborted.')
                alignment = ''
            sep = line.rfind('_')
            OG = line[1:sep]
            seqid = line[sep + 1:]
            species = seqid[:4]
            species_seqids_dict[species] = species_seqids_dict.get(species, []) + [seqid]
            OG_dict[OG] = OG_dict.get(OG, []) + [seqid]
            if species not in supergroups:
                sys.stderr.write(f"Warning: {species} not recognised.\n")
        else:
            alignment += line
    if aln_length is None:
        sys.exit('Error: no sequences found in the alignment file! Analysis aborted.')
    elif len(alignment) != aln_length:
        sys.exit('Error: not all aligned sequences have the same length! Analysis aborted.')
    aligned_proteins[seqid] = alignment
    return species_seqids_dict, OG_dict, aligned_proteins

def get_coordinates(species_seqids_dict, euk_path = '/home/julian/julian2/snel-clan-genomes/eukarya_new'):
    seqid_coordinates = {}
    for species, seqids in species_seqids_dict.items():
        if species in ('BBRI', 'ESIL'):
            sys.stderr.write(f"Warning: no CDS coordinates can be retrieved from the GFF file for {species}.\n")
            continue
        if os.path.isfile(f'{euk_path}/data_set/gff_files/{species}.gff3.gz'):
            gff_file_name = f'{euk_path}/data_set/gff_files/{species}.gff3.gz'
        elif os.path.isfile(f'{euk_path}/data_set/gff_files/{species}.gff.gz'):
            gff_file_name = f'{euk_path}/data_set/gff_files/{species}.gff.gz'
        else:
            sys.stderr.write(f'Warning: no GFF file found for {species}.\n')
            continue
        original_ids = {}
        sys.stderr.write(f'Reading metadata file for {species}...\n')
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
        sys.stderr.write(f'Reading GFF file for {species}...\n')
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
                        sys.stderr.write(f'Error: elements in this line in the GFF file for {species} not recognised:\n{line}\n')
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
                        new_entry = [int(fields[3]), int(fields[4]), fields[6], fields[7]]
                        try:
                            seqid_coordinates[eukarya_id].append(new_entry)
                        except KeyError:
                            seqid_coordinates[eukarya_id] = [new_entry]
        for seqid in seqids:
            if seqid not in seqid_coordinates:
                sys.stderr.write(f'Warning: {seqid} not detected in GFF file.\n')
        sys.stderr.write('Done!\n')
    return seqid_coordinates

def map_introns_aa(euk_cds_dict, lengths):
    location_introns = {}
    for seqid, CDSs in euk_cds_dict.items():
        # Check if all coding parts of one gene lie in the same direction
        directions = set([cds[-2] for cds in CDSs])
        if len(directions) != 1:
            sys.stderr.write(f"Warning: the CDSs of {seqid} are not in the same direction. Excluded from analysis.\n")
        location_introns[seqid] = []
        length_without_introns = 0
        #location_introns = []
        startcds = CDSs[0][0]
        stopcds = CDSs[-1][1] #-1 needed for counting in python
        # Some of the minus directed CDSs needed to be reversed.
        if directions == {"-"}:
            if startcds < stopcds:
                CDSs.reverse()
        try:
            start_phase = int(CDSs[0][3])
        except ValueError:
            sys.stderr.write(f'Warning: no start phase detected for {seqid}. Excluded from analysis.\n')
        if start_phase != 0:
            if directions == {'-'}:
                CDSs[0][1] += 3 - start_phase
            else:
                CDSs[0][0] -= 3 - start_phase
        # Loop through CDSs, excluding the last to prevent the end of protein being seen as intron, and add to list.
        for intron_information in CDSs[:-1]:
            length_without_introns += intron_information[1] - intron_information[0] + 1
            phase = length_without_introns % 3
            location_intron = int((length_without_introns - phase)/3 + 1)
            location_introns[seqid].append([phase, location_intron])
        length_without_introns += CDSs[-1][1] - CDSs[-1][0] + 1
        # everything is devided by three, from length of mRNA (nucleotides) to polypeptide length (amino acid).
        if length_without_introns%3 != 0:
            sys.stderr.write(f"Warning: {seqid} seems to be incorrectly annotated.\n") #A small check if the genes are correctly annotated
        if not (length_without_introns // 3 == lengths[seqid] or length_without_introns // 3 == lengths[seqid] + 1):
            sys.stderr.write(f'Warning: different lengths for {seqid}: {length_without_introns // 3} (gff) | {lengths[seqid]} (aln). Excluded from analysis.\n')
            del location_introns[seqid]
    return location_introns

def count_groups(seqids, supergroups, unique = True):
    group_counts = {g : 0 for g in set(supergroups.values())}
    species = [seqid[:4] for seqid in seqids]
    if unique:
        species = set(species)
    for spec in species:
        group_counts[supergroups[spec]] += 1
    return group_counts

def map_introns_aln(location_introns, aligned_proteins, seqid_OG):
    aln_intron_positions = {}
    for seqid, introns in location_introns.items():
        if len(introns) == 0:
            continue
        intron_count = 0 #regulates that not every intron is checked every time.
        location = 0
        og = seqid_OG[seqid]
        if og not in aln_intron_positions:
            aln_intron_positions[og] = {}
        og_positions = aln_intron_positions[og]  # Retrieve information present
        for position, character in enumerate(aligned_proteins[seqid]):
            if character == "-":
                continue
            location += 1
            # Multiple introns could be located on this amino acid
            while introns[intron_count][1] == location:
                phase = introns[intron_count][0]
                aln_intron_position = position + 1
                if aln_intron_position not in og_positions:
                    og_positions[aln_intron_position] = [[], [], []]
                og_positions[aln_intron_position][phase].append(seqid)
                intron_count += 1
                if intron_count == len(introns):
                    break
            else:
                continue
            break
        else:
            sys.stderr.write(f'Not all intron positions for {seqid} found in the alignment.\n')
    return aln_intron_positions

def cluster_neighbouring_positions(aln_intron_positions, OG_dict, nt_shifts, aa_shifts):
    introns_shifts = {}
    for OG, aa_positions in aln_intron_positions.items():
        introns_shifts[OG] = {}
        for position, phases in aa_positions.items():
            introns_shifts[OG][position] = {}
            for phase, seqs in enumerate(phases):
                if len(seqs) == 0:
                    continue
                introns_shifts[OG][position][phase] = []
                for k in range(int(position-aa_shifts), int(position+aa_shifts) + 1, 1):  #k is just a variable looping through the list of numbers
                    if k in aln_intron_positions[OG]:
                        introns_shifts[OG][position][phase].extend(aln_intron_positions[OG][k])
                    else:
                        introns_shifts[OG][position][phase].extend([[],[],[]])
                introns_shifts[OG][position][phase] = introns_shifts[OG][position][phase][int(3*aa_shifts+phase-nt_shifts):int(3*aa_shifts+phase+nt_shifts+1)]
                for j in range(0, len(introns_shifts[OG][position][phase]), 1):
                    if len(introns_shifts[OG][position][phase][j]) > len(seqs):
                        #Assumption is that in the range defined there can only be 1 LECA intron and that if shift <3
                        #there can only be one dominant intron on 1 aa position.
                        del introns_shifts[OG][position][phase]
                        #Selection: If near this position is a location with more introns, this location is deleted.
                        break
            if len(introns_shifts[OG][position]) == 0:
                del introns_shifts[OG][position]
    return introns_shifts

def get_leca_introns(og_introns, supergroups, OG_group_species_count, OG_group_seq_count, threshold_percentage_genes, threshold_percentage_species, threshold_species, outdir):
    leca_introns = {}
    leca_count = {}
    for OG, positions in og_introns.items():
        leca_count[OG] = 0
        og_file = open(f'{outdir}/{OG}_introns.tsv', 'w')
        group_order = sorted(set(supergroups.values()))
        og_file.write('Position\tPhase\tLECA\tTotal species (%)\t')
        og_file.write('\t'.join([group + ' species' for group in group_order]) + '\t')
        og_file.write('Total sequences (%)\t')
        og_file.write('\t'.join([group + ' sequences' for group in group_order]) + '\n')
        spec_total = sum(OG_group_species_count[OG].values())
        seq_total = sum(OG_group_seq_count[OG].values())
        leca_introns[OG] = {}
        for position, phases in positions.items():
            for phase, seqids in phases.items():
                seqids = [seqid for seqids in introns_shifts[OG][position][phase] for seqid in seqids]
                group_species_counts = count_groups(seqids, supergroups, unique = True)
                opimoda = group_species_counts['Amoebozoa'] + group_species_counts['Obazoa']
                diphoda = sum(group_species_counts.values()) - opimoda
                spec_cov = sum(group_species_counts.values()) / spec_total * 100
                group_seq_counts = count_groups(seqids, supergroups, unique = False)
                seq_cov = sum(group_seq_counts.values()) / seq_total * 100
                leca = 'No'
                if opimoda >= threshold_species and diphoda >= threshold_species:
                    if spec_cov >= threshold_percentage_species and seq_cov >= threshold_percentage_genes:
                        leca = 'Yes'
                        leca_count[OG] += 1
                        if position in leca_introns[OG]:
                            leca_introns[OG][position].append(phase)
                        else:
                            leca_introns[OG][position] = [phase]
                print(position, phase, leca, round(spec_cov, 1), sep = '\t', end = '', file = og_file)
                for group, total in OG_group_species_count[OG].items():
                    if total != 0:
                        og_file.write(f'\t{round(group_species_counts[group] / total * 100, 1)}')
                    else:
                        og_file.write('\tNA')
                for group, total in OG_group_seq_count[OG].items():
                    if total != 0:
                        og_file.write(f'\t{round(group_seq_counts[group] / total * 100, 1)}')
                    else:
                        og_file.write('\tNA')
                og_file.write('\n')
        og_file.close()
    return leca_introns, leca_count

def get_shared_leca_introns(leca_introns, nt_shifts, aa_shifts):
    OGs = leca_introns.keys()
    number_of_shared_introns = {}
    for i in range(len(OGs) - 1):
        OG1 = tuple(OGs)[i]
        number_of_shared_introns[OG1] = {}
        for j in range(i + 1, len(OGs)):
            OG2 = tuple(OGs)[j]
            number_of_shared_introns[OG1][OG2] = 0
            OG2_introns_copy = copy.deepcopy(leca_introns[OG2])
            for position, phases in leca_introns[OG1].items():
                for phase1 in phases:
                    for possible_shift in range(position - aa_shifts, position + aa_shifts + 1):
                        if possible_shift in OG2_introns_copy:
                            for phase2 in OG2_introns_copy[possible_shift]:
                                nt = possible_shift*3 + phase2
                                if nt in range(position*3 + phase1 - nt_shifts, position*3 + phase1 + nt_shifts + 1):
                                    number_of_shared_introns[OG1][OG2] += 1
                                    sys.stderr.write(f'Found a shared intron between {OG1} ({position}.{phase1}) and {OG2} ({possible_shift}.{phase2}).\n')
                                    OG2_introns_copy[possible_shift].remove(phase2)
                                    if len(OG2_introns_copy[possible_shift]) == 0:
                                        del OG2_introns_copy[possible_shift]
                                    break
                            else:
                                continue
                            break
    return number_of_shared_introns

# Parse arguments
parser = argparse.ArgumentParser(description = "This script maps intron positions onto a protein alignment.")
parser.add_argument("alignment", help = 'protein alignment in fasta format')
parser.add_argument('-o', metavar = 'outdir', help = 'directory for output files (default: current)')
parser.add_argument('-e', metavar = 'eukarya_path', help = 'directory containing the Eukarya database (default: ~julian/julian2/snel-clan-genomes/eukarya)')
parser.add_argument('-i', help = 'infer LECA introns and introns predating duplications (default: off)', action = 'store_true')
parser.add_argument('-s', metavar = 'shifts', help = 'number of nucleotide shifts allowed for an intron (default: 0)', type = int, default = 0)
parser.add_argument("-p", metavar = 'gene%', help = "percentage of genes in which an intron at least has to occur to call it a LECA intron (default: 7.5)", default = 7.5)
parser.add_argument("-t", metavar = 'species%', help = "percentage of species in which an intron at least has to occur to call it a LECA intron (default: 15)", default = 15)
parser.add_argument("-n", metavar = 'species', help = "minimum number of Opimoda and Diphoda species that should have an intron to call it a LECA intron (default: 2)", default = 2)
#parser.add_argument('-d', metavar = 'domains', help = 'Pfam domains instead of full-length genes', action = 'store_true')
args = parser.parse_args()

# To add: sequence ID --> paralogue table

# Parse arguments
prefix = os.path.basename(args.alignment).split('.')[0]
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
    sys.stderr.write(f"Warning: the number of shifts you want to tolerate ({shifts}) seems very high\n.")

# Set thresholds for considering an intron a LECA intron
threshold_percentage_genes = args.p # Intron has to be present in at least this many genes
threshold_percentage_species = args.t # Intron has to be present in at least this many species
threshold_species = args.n # Intron has to be present in at least X Opimoda (~unikont) and X Diphoda (~bikont) species

info = f'Intron maper v1.0\nDeveloped by Sjoerd Gremmen, Michelle Schinkel and Julian Vosseberg.\nTime: {time.asctime()}\nCommand: {" ".join(sys.argv)}\n\n'
if inference:
    info += f'For LECA inference:\n-Shifts: {shifts}\n-Gene%: {threshold_percentage_genes}\n-Species%: {threshold_percentage_species}\n-Opimoda/Diphoda: {threshold_species}\n\n'
sys.stderr.write(info)
log = open(f'{output_path}/{prefix}.log', 'w')
log.write(info)

groups = set(supergroups.values())

# Step 1: Parse input files

# Parse alignment
sys.stderr.write(f'Reading alignment file {args.alignment}...\n')
with open(args.alignment) as fasta_file:
    species_seqids_dict, OG_dict, aligned_proteins = parse_alignment(fasta_file)
lengths = {seqid : len(aln.replace('-', '')) for seqid, aln in aligned_proteins.items()}
info = f'Alignment has {len(aligned_proteins)} sequences from {len(OG_dict)} OGs with {len(tuple(aligned_proteins.values())[0])} columns.\n'
info += "\n".join([OG + ": " + str(len(seqids)) + " sequences" for OG, seqids in OG_dict.items()]) + '\n\n'
sys.stderr.write(info)
log.write(info)

# Get CDS coordinates
sys.stderr.write('Obtaining CDS coordinates...\n')
euk_cds_dict = get_coordinates(species_seqids_dict)
info = f'Intron location information for {len(euk_cds_dict)} / {len(aligned_proteins)} sequences.\n\n'
sys.stderr.write(info)
log.write(info)

# Step 2: Map intron positions onto the proteins
sys.stderr.write('Mapping intron positions onto the proteins...\n')
location_introns = map_introns_aa(euk_cds_dict, lengths)
info = f'Mapping successful for {len(location_introns)} / {len(euk_cds_dict)} sequences.\n\n'
sys.stderr.write(info)
log.write(info)
OG_seq_info = {}
for OG, seqs in OG_dict.items():
    OG_seq_info[OG] = []
    for seq in seqs:
        if seq in location_introns:
            OG_seq_info[OG].append(seq)
OG_group_seq_count = {OG : count_groups(seqs, supergroups, unique = False) for OG, seqs in OG_seq_info.items()}
OG_group_species_count = {OG : count_groups(seqs, supergroups, unique = True) for OG, seqs in OG_seq_info.items()}

# Step 3: Create sequence features file for viewing intron positions in Jalview
sys.stderr.write('Creating Jalview sequence features file...\n')
seqid_OG = {seqid : og for og, seqids in OG_dict.items() for seqid in seqids}
with open(f'{output_path}/{prefix}_jalview.sff', 'w') as features_file:
    features_file.write('phase0\tgreen\nphase1\tblue\nphase2\tmagenta\nNA\tblack\n')
    for seqid in aligned_proteins:
        try:
            introns = location_introns[seqid]
            for phase, position in introns:
                features_file.write(f'Intron position\t{seqid_OG[seqid]}_{seqid}\t-1\t{position}\t{position}\tphase{phase}\t\n')
        except KeyError:
            features_file.write(f'Intron position\t{seqid_OG[seqid]}_{seqid}\t-1\t1\t{lengths[seqid]}\tNA\t\n')
sys.stderr.write('Done!\n\n')

# Step 4: Map intron positions onto the alignment
sys.stderr.write('Mapping intron positions onto the alignment...\n')
aln_intron_positions = map_introns_aln(location_introns, aligned_proteins, seqid_OG)
with open(f'{output_path}/{prefix}_introns.tsv', 'w') as all_introns:
    all_introns.write('Positon\tPhase\tOG\tSequence IDs\n')
    for OG, positions in aln_intron_positions.items():
        for position, phases in positions.items():
            for phase, seqids in enumerate(phases):
                if len(seqids) == 0:
                    continue
                print(position, phase, OG, ','.join(seqids), sep = '\t', file = all_introns)
sys.stderr.write(f'Written to {prefix}_introns.tsv!\n\n')

# Step 5: infer LECA introns
if inference:
    sys.stderr.write('Inferring introns in LECA and reporting LECA introns shared between OGs...\n')
    # Calculate shifts and cluster neighbouring positions
    if shifts <= 3:
        aa_shifts = 1
    else:
        if shifts%3 == 0:
            aa_shifts = shifts/3
        else:
            aa_shifts =(shifts - shifts%3)/3 + 1 # 0-3 --> aa=1, 4-6 --> aa=2, 7-9 --> aa=3, ...
    introns_shifts = cluster_neighbouring_positions(aln_intron_positions, OG_dict, shifts, aa_shifts)

    # Infer intron positions/clusters present in LECA and make an analysis file per OG
    leca_introns, leca_count = get_leca_introns(introns_shifts, supergroups, OG_group_species_count, OG_group_seq_count, threshold_percentage_genes, threshold_percentage_species, threshold_species, output_path)
    info = ' LECA intron(s) inferred.\n'
    total = 0
    for OG, leca_intron_count in leca_count.items():
        info += f'{OG}: {leca_intron_count} LECA intron(s)\n'
        total += leca_intron_count
    info = str(total) + info
    sys.stderr.write(info)
    log.write(info)

    # Obtain shared LECA introns
    number_of_shared_introns = get_shared_leca_introns(leca_introns, shifts, aa_shifts)
    info = f' / {total} LECA intron(s) shared between OGs.\n'
    total = 0
    for OG1, OG2s in number_of_shared_introns.items():
        for OG2, count in OG2s.items():
            if count > 0:
                info += f'{OG1} - {OG2}: {count} shared LECA intron(s)\n'
                total += count
    info = str(total) + info
    sys.stderr.write(info)
    log.write(info)

    # Write output table
    with open(output_path + '/table.txt', 'w') as table:
        table.write('\t' + '\t'.join(OG_dict.keys()) + '\n')
        OGs = list(OG_dict.keys())
        for i, OG1 in enumerate(OGs):
            table.write(OG1)
            for j in range(i + 1):
                OG2 = OGs[j]
                if OG1 == OG2:
                    table.write(f'\t{leca_count[OG1]}')
                else:
                    if leca_count[OG1] == 0:
                        table.write('\t0')
                        continue
                    try:
                        number_shared = number_of_shared_introns[OG1][OG2]
                    except KeyError:
                        number_shared = number_of_shared_introns[OG2].get(OG1, 0)
                    table.write(f'\t{number_shared}')
            table.write('\n')
log.close()
