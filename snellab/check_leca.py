#!/usr/bin/env python3
from eukarya import *
import sys
from numpy import median

def assign_all_seqs(bbh_seqs, all_seqs):
    """Assigns all original sequences to their best representing sequence based on the BLAST results and keeps track of human sequences"""
    other_seqs = set(all_seqs) - set(bbh_seqs)
    human_represent = {}
    for seq in bbh_seqs:
        if 'HSAP' in seq:
            human_represent[seq] = seq # Human seq represented by itself
    blast = open('/home/julian/julian2/pfam_hmm/improved_pipeline/euk_fasta/blast_output/' + prefix + '/' + prefix + '_blastp.txt')
    best_blast_hits = {}
    for line in blast:
        line = line.rstrip()
        fields = line.split('\t')
        query = fields[0]
        if query in bbh_seqs:
            continue
        hit = fields[1]
        if hit not in bbh_seqs:
            continue
        bit_score = float(fields[11])
        if query not in best_blast_hits:
            best_blast_hits[query] = (hit, bit_score)
        else:
            best_score = best_blast_hits[query][1]
            if bit_score > best_score:
                best_blast_hits[query] = (hit, bit_score)
    blast.close()
    representing = {}
    for query in best_blast_hits:
        best_hit = best_blast_hits[query][0]
        if best_hit in representing:
            representing[best_hit].append(query)
        else:
            representing[best_hit] = [query]
        if 'HSAP' in query:
            human_represent[query] = best_hit
    return representing, human_represent

def infer_coverage_redundancy(node, representing):
    number_species = {'Obazoa':123, 'Amoebozoa':6, 'RASH':43, 'Archaeplastida':24, 'Excavata':13}
    all_species = []
    species_counts = {}
    supergroup_counts = {'Obazoa':0, 'Amoebozoa':0, 'RASH':0, 'Archaeplastida':0, 'Excavata':0}
    group_species_counts = {'Obazoa':[], 'Amoebozoa':[], 'RASH':[], 'Archaeplastida':[], 'Excavata':[]}
    for leaf in node: # node can also be a list of leaves instead
        name = leaf
        sp = name[0:4]
        species = [seq[0:4] for seq in representing.get(name,'')]
        species.append(sp)
        all_species.extend(species)
    for spec in all_species:
        if spec in species_counts:
            species_counts[spec] += 1
        else:
            supergroup = supergroups5[spec]
            supergroup_counts[supergroup] += 1
            species_counts[spec] = 1
    coverages = []
    for group, counts in supergroup_counts.items():
        coverages.append(counts / number_species[group])
    coverage_av = sum(coverages) / 5
    for spec, counts in species_counts.items():
        group = supergroups5[spec]
        group_species_counts[group].append(counts)
    redundancy_med = median([median(counts) for counts in group_species_counts.values() if len(counts) > 0])
    return coverage_av, redundancy_med, list(species_counts.keys())

def get_human_representing(seqs, human_seqs, human_represent, human_seq_info):
    remaining_human_seqs = human_seqs[:]
    repr_human = {}
    for human in human_seqs:
        if human in human_represent:
            if human_represent[human] in seqs:
                if '_' in human:
                    human_seq = human[:human.find('_')]
                else:
                    human_seq = human
                repr_human[human] = human_seq_info[human_seq][0]
                remaining_human_seqs.remove(human)
    return remaining_human_seqs, repr_human

# ------

coverage_criterion = 0.15
prefix = sys.argv[1]
supergroups2, supergroups5 = get_supergroups()
bbh_seqs = sys.argv[2:]
all_seqs_file = open('/home/julian/julian2/pfam_hmm/improved_pipeline/euk_fasta/' + prefix + '_seqids.list')
all_seqs = [line.rstrip() for line in all_seqs_file]
all_seqs_file.close()
representing, human_represent = assign_all_seqs(bbh_seqs, all_seqs)
human_seq_info = seq_info(['HSAP'])['HSAP']
human_seqs = [seq for seq in all_seqs if 'HSAP' in seq]
coverage, copy_no, species = infer_coverage_redundancy(bbh_seqs, representing)
if coverage <= coverage_criterion:
    sys.exit('No LECA in this Pfam. Coverage is %f.' % coverage)
else:
    lecas_out = open(prefix + '_lecas.tsv', 'w')
    print('Pfam\tFECA\tAncestry\tLECA\tSupport\tCoverage\tCopy number\tHuman seqs\tHuman name\tSeqs', file = lecas_out)
    if len(human_represent) > 0:
        human_seqs, leca_human = get_human_representing(bbh_seqs, human_seqs, human_represent, human_seq_info)
        print(prefix, 'NA', 'Eukaryotic', 'OG1.1', 'NA', coverage, copy_no, ','.join(list(leca_human.keys())), ','.join(list(leca_human.values())), ','.join(bbh_seqs), sep = '\t', file = lecas_out)
    else:
        print(prefix, 'NA', 'Eukaryotic', 'OG1.1', 'NA', coverage, copy_no, '', '', ','.join(bbh_seqs), sep = '\t', file = lecas_out)
    lecas_out.close()
    print('Pfam\tFECAs (normal)\tFECAs (strict)\tFECAs (after merging)\tLECAs\tUnknowns\tNon-FECAs')
    print(prefix, 'NA', 'NA', 'NA', 1, 0, 'NA', sep = '\t')
