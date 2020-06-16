#!/usr/bin/env python3

import sys
import os
from subprocess import PIPE, Popen, call

pfam = sys.argv[1]
groups = sys.argv[2]
if len(sys.argv) > 3:
    outdir = sys.argv[3]
else:
    outdir = '.'

if groups == "2":
    path_to_files = '/home/julian/julian2/timing_dupl/4_phylogeny_pfam_pipeline2/1_bbhs_2_groups/4_ete_analysis/d20_l15/separate/' + pfam + '/'
elif groups == "5":
    path_to_files = '/home/julian/julian2/timing_dupl/4_phylogeny_pfam_pipeline2/2_bbhs_5_groups/4_ete_analysis/d20_l15/separate/' + pfam + '/'
elif groups == 'r2':
    path_to_files = '/home/julian/julian2/timing_dupl/6_nee_revision/1_kclust_incl_asgard/5_ete_analysis/bbhs_2_groups/' + pfam + '/'
    if not os.path.exists(path_to_files):
        path_to_files = '/home/julian/julian2/timing_dupl/4_phylogeny_pfam_pipeline2/1_bbhs_2_groups/4_ete_analysis/d20_l15/separate/' + pfam + '/'
else:
    sys.exit("Error: groups not recognised!")
fecas = {}

lecas_file = open(path_to_files + pfam + '_lecas.tsv')
for line in lecas_file:
    if line.startswith('Pfam'):
        continue
    line = line.rstrip()
    fields = line.split('\t')
    feca = fields[1]
    seqs = fields[9]
    try:
        fecas[feca] += ',' + seqs
    except KeyError:
        fecas[feca] = seqs
lecas_file.close()

unknowns_file = open(path_to_files + pfam + '_unknowns.tsv')
for line in unknowns_file:
    if line.startswith('Pfam'):
        continue
    line = line.rstrip()
    fields = line.split('\t')
    feca = fields[1]
    seqs = fields[7]
    fecas[feca] += ',' + seqs
unknowns_file.close()

mfecas = {}
fecas_file = open(path_to_files + pfam + '_fecas_strict.tsv')
for line in fecas_file:
    if line.startswith('Pfam'):
        continue
    line = line.rstrip()
    fields = line.split('\t')
    mfeca = fields[1]
    if '_' in mfeca:
        content = fields[5].split('+')
        mfecas[mfeca] = content
    else:
        mfecas[mfeca] = [mfeca]
fecas_file.close()

mfecas_seqs = {}
for mfeca, content in mfecas.items():
    for feca in content:
        if 'Non-FECA' in feca:
            seqs = Popen(f'grep -w {feca} {path_to_files}{pfam}_non_feca.tsv | cut -f 8', stdout = PIPE, shell = True).communicate()[0].decode().rstrip()
            try:
                mfecas_seqs[mfeca] += ',' + seqs
            except KeyError:
                mfecas_seqs[mfeca] = seqs
        else:
            try:
                mfecas_seqs[mfeca] += ',' + fecas[feca]
            except KeyError:
                mfecas_seqs[mfeca] = fecas[feca]

if groups == "2" or groups == 'r2':
    fasta = '/home/julian/julian2/pfam_hmm/improved_pipeline/combined/2_groups/' + pfam + '_2_groups.fa'
else:
    fasta = '/home/julian/julian2/pfam_hmm/improved_pipeline/combined/5_groups/' + pfam + '_5_groups.fa'
for mfeca, seqs in mfecas_seqs.items():
    call(f'selectSequences.pl -t -e -i {fasta} -q {seqs} > {outdir}/{pfam}_{mfeca}.fa', shell = True)
