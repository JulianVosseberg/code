#!/usr/bin/env python3

import sys

def usage():
    sys.exit("""\n\tUsage: leca_pfams_to_leca_genes.py <files>
    
This script identifies the overlap between LECA OGs.""")

if len(sys.argv) == 1:
    usage()
else:
    files = sys.argv[1:]

overlap_cutoff = 10
consistency_cutoff = 0.2
cutoff = 'consistency'

#lecas_all_seqs_file = open('/home/julian/julian2/timing_dupl/4_phylogeny_pfam_pipeline2/1_bbhs_2_groups/3_iqtree_LG4X/test/PF00004_lecas_all_seqs.tsv')
lecas = {}
for filename in files:
    lecas_all_seqs_file = open(filename)
    for line in lecas_all_seqs_file:
        if line.startswith('Pf'):
            continue
        line = line.rstrip()
        fields = line.split('\t')
        leca = fields[0] + '_' + fields[1]
        seqs = fields[2].split(',') + fields[3].split(',')
        for i,seq in enumerate(seqs):
            if '_' in seq:
                seqs[i] = seq[:seq.find('_')]
        lecas[leca] = set(seqs)
    lecas_all_seqs_file.close()

lecas_list = list(lecas.keys())
for i in range(len(lecas_list) - 1):
    for j in range(i + 1, len(lecas_list)):
        leca1_seqs = lecas[lecas_list[i]]
        leca2_seqs = lecas[lecas_list[j]]
        overlap = len(leca1_seqs & leca2_seqs)
        consistency = overlap / len(leca1_seqs | leca2_seqs)
        if cutoff == 'consistency':
            if consistency < consistency_cutoff:
                continue
        elif cutoff == 'overlap':    
            if overlap < 10:
                continue
        else:
            sys.exit('Cutoff not recognised')
        print(lecas_list[i], lecas_list[j], overlap, consistency)
