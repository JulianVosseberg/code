#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 14:22:51 2024

@author: julian
"""

import sys

taxon_seqs = {}
with open(sys.argv[1]) as aln:
    seq = ''
    for line in aln:
        line = line.rstrip()
        if line.startswith('>'):
            if seq != '':
                taxon_seqs[taxon] = seq
            taxon = line[1:]
            seq = ''
        else:
            seq += line.replace('-', '').replace('X', '').replace('*', '')
        taxon_seqs[taxon] = seq

total = ''.join(taxon_seqs.values())
if set(total) == set('ACGT'):
    alphabet = 'ACGT'
else:
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
freqs = {letter: total.count(letter) / len(total) for letter in alphabet}
chi2_total = 0
for taxon, seq in taxon_seqs.items():
    chi2 = 0
    for letter in alphabet:
        expected = freqs[letter]*len(seq)
        chi2 += (seq.count(letter)-expected)**2 / expected
    print(taxon, round(chi2, 2))
    chi2_total += chi2
print('\nTotal', round(chi2_total, 2))