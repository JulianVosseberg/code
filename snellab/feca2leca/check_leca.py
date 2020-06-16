#!/usr/bin/env python3
import sys
sys.path.append('/home/julian/julian2/scripts/github/snellab/')
from eukarya import *
from scrollsaw import *
import argparse

# ------

parser = argparse.ArgumentParser(description = "This script checks if a set of sequences fulfills the LECA criterion by using the corresponding BLAST hits.")
parser.add_argument("pfam", metavar = "Pfam", help = "the Pfam for which you want to check if it was in LECA")
parser.add_argument("sequence_ids", metavar = "sequence_ID", nargs = '+', help = 'the sequences from the BBHs')
parser.add_argument("-l", metavar = "0.xx", help = "coverage threshold for LECA calling (DEFAULT: 0.15)", type = float, default = 0.15)
parser.add_argument("-r", metavar = "root", help = "position of eukaryotic root (DEFAULT: Opimoda-Diphoda)", default = "OD")
parser.add_argument("-s", metavar = "supergroups", help = "supergroups definition used", type = int, choices = (4, 5, 6), default = 5)
parser.add_argument('-o', metavar = 'outdir', help = 'directory for output files (DEFAULT: current)', default = '.')
args = parser.parse_args()

prefix = args.pfam
bbh_seqs = args.sequence_ids
coverage_criterion = args.l
outdir = args.o
if args.s == 5:
    supergroups = supergroups5
elif args.s == 4:
    supergroups = supergroups4
else:
    supergroups = supergroups6
root_daughters = get_root_daughters(args.r, supergroups4)
root_groups = set([root_daughters[seq[:4]] for seq in bbh_seqs])
if len(root_groups) != 2:
    sys.exit(f'No LECA in this Pfam. Only {tuple(root_groups)[0]} BBH sequences.')

with open(f'/home/julian/julian2/pfam_hmm/improved_pipeline/euk_fasta/{prefix}_seqids.list') as all_seqs_file:
    all_seqs = [line.rstrip() for line in all_seqs_file]
representing, human_represent = assign_all_seqs(bbh_seqs, all_seqs, euk_only = True, prefix = prefix)
human_seq_info = seq_info(['HSAP'])['HSAP']
human_seqs = [seq for seq in all_seqs if 'HSAP' in seq]
coverage, copy_no, species = infer_coverage_redundancy(bbh_seqs, representing, supergroups, tree = False)
if coverage <= coverage_criterion:
    sys.exit(f'No LECA in this Pfam. Coverage is {coverage}.')
else:
    with open(f'{outdir}/{prefix}_lecas.tsv', 'w') as lecas_out:
        lecas_out.write('Pfam\tFECA\tAncestry\tLECA\tSupport\tCoverage\tCopy number\tHuman seqs\tHuman name\tSeqs\n')
        if len(human_represent) > 0:
            human_seqs, leca_human = get_human_representing(bbh_seqs, human_seqs, human_represent, human_seq_info)
            print(prefix, 'NA', 'Eukaryotic', 'OG1.1', 'NA', coverage, copy_no, ','.join(leca_human.keys()), ','.join(leca_human.values()), ','.join(bbh_seqs), sep = '\t', file = lecas_out)
        else:
            print(prefix, 'NA', 'Eukaryotic', 'OG1.1', 'NA', coverage, copy_no, '', '', ','.join(bbh_seqs), sep = '\t', file = lecas_out)
    print('Pfam\tFECAs (normal)\tFECAs (strict)\tFECAs (after merging)\tLECAs\tUnknowns\tNon-FECAs')
    print(prefix, 'NA', 'NA', 'NA', 1, 0, 'NA', sep = '\t')
