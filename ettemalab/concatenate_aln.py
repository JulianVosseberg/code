#!/usr/bin/env python3

import glob
import sys
import os
import argparse

def add_to_supermatrix(alignments, marker, taxon, seq, aln_file):
    if taxon not in alignments:
        alignments[taxon] = {}
    if marker in alignments[taxon]:
        sys.exit(f'Error: multiple hits for {taxon} in {aln_file}.')
    alignments[taxon][marker] = seq

# Parse arguments
parser = argparse.ArgumentParser(description = "This script concatenates single alignments into one alignment (supermatrix), which is printed to stdout.")
parser.add_argument("alignments", help = "single alignments in FASTA format", nargs = "*")
parser.add_argument('-p', metavar = "prefix", help = "prefix delimiter")
parser.add_argument('-s', metavar = "suffix", help = 'suffix delimiter')
args = parser.parse_args()


marker_len = {}
alignments = {}
for marker, aln_file in enumerate(args.alignments):
    with open(aln_file) as aln:
        seq = ''
        for line in aln:
            line = line.rstrip()
            if line.startswith('>'):
                if seq != '':
                    add_to_supermatrix(alignments, marker, taxon, seq, aln_file)
                taxon = line[1:]
                if args.p:
                    posit = taxon.find(args.p)
                    if posit != -1:
                        taxon = taxon[posit+1:]
                if args.s:
                    posit = taxon.find(args.s)
                    if posit != -1:
                        taxon = taxon[:posit]
                seq = ''
            else:
                seq += line
        add_to_supermatrix(alignments, marker, taxon, seq, aln_file)
        marker_len[marker] = len(seq)

supermatrix = {}
for taxon, marker_alns in alignments.items():
    superaln = ''
    for marker in marker_len:
        if marker not in marker_alns:
            superaln += marker_len[marker] * '-'
        else:
            superaln += marker_alns[marker]
    supermatrix[taxon] = superaln
 
for taxon, superaln in supermatrix.items():
    sys.stdout.write(f'>{taxon}\n')
    for i in range(0, len(superaln), 70):
        sys.stdout.write(f'{superaln[i:i+70]}\n')
