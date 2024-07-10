#!/usr/bin/env python3

import argparse
import sys
from Bio import AlignIO
import numpy as np

def calculate_chi2_scores(alignment, report = True):
    taxon_seqs = {}
    for record in alignment:
        seq = record.seq.replace('-', '').replace('X', '').replace('*', '')
        taxon_seqs[record.id] = str(seq)
    total = ''.join(taxon_seqs.values())
    if set(total) == set('ACGT'):
        alphabet = 'ACGT'
        critical_value = 7.814728
    else:
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        critical_value = 30.143527
    freqs = {letter: total.count(letter) / len(total) for letter in alphabet}
    chi2_total = 0
    failed = 0
    for taxon, seq in taxon_seqs.items():
        chi2 = 0
        for letter in alphabet:
            expected = freqs[letter]*len(seq)
            chi2 += (seq.count(letter)-expected)**2 / expected
        if chi2 > critical_value:
            sign = '(*)'
            failed += 1
        else:
            sign = ''
        if report:
            print(taxon, round(chi2, 2), sign)
        chi2_total += chi2
    print('\nTotal', round(chi2_total, 2), f'({failed} sequences failed the composition chi2 test (p < 0.05))')
    return chi2_total, alphabet

def calculate_deltachi2_positions_old(alignment, chi2_total, alphabet):
    length = alignment.get_alignment_length()
    position_deltachi2 = {}
    for i in range(length):
        sys.stderr.write(f"\rProcessing column {i+1:>6}...")
        alignment_wo_i = alignment[:,:i] + alignment[:,i+1:]
        taxon_seqs = {}
        for record in alignment_wo_i:
            seq = record.seq.replace('-', '').replace('X', '').replace('*', '')
            taxon_seqs[record.id] = str(seq)
        total = ''.join(taxon_seqs.values())
        freqs = {letter: total.count(letter) / len(total) for letter in alphabet}
        chi2 = 0
        for seq in taxon_seqs.values():
            for letter in alphabet:
                expected = freqs[letter]*len(seq)
                chi2 += (seq.count(letter)-expected)**2 / expected
        position_deltachi2[i] = chi2 - chi2_total
    return position_deltachi2

def calculate_deltachi2_positions(alignment, freqs, total_length, taxon_seqs, chi2_total, alphabet):
    length = alignment.get_alignment_length()
    position_deltachi2 = {}
    for i in range(length):
        sys.stderr.write(f"\rProcessing column {i+1:>6}...")
        column = alignment[:,i]
        column_wo_gaps = column.replace('-', '').replace('X', '').replace('*', '')
        freqs_wo_i = {letter: (freqs[letter] - column_wo_gaps.count(letter)) / (total_length - len(column_wo_gaps)) for letter in alphabet}
        chi2 = 0
        for seq in taxon_seqs:
            seq_length = len(seq.replace('-', '').replace('X', '').replace('*', ''))
            if seq[i] not in ('-', 'X', '*'):
                seq_length -= 1
            for letter in alphabet:
                expected = freqs_wo_i[letter]*seq_length
                if seq[i] == letter:
                    chi2 += (seq.count(letter)-1-expected)**2 / expected
                else:
                    chi2 += (seq.count(letter)-expected)**2 / expected
        position_deltachi2[i] = chi2 - chi2_total
    return position_deltachi2

def parse_blocks(alignment, blocksfile):
    length = alignment.get_alignment_length()
    blocks = []
    i = 0
    line = blocksfile.readline().rstrip()
    fields = line.split()
    if len(fields) != 3:
        sys.exit('Error: blocks file does not contain three fields (<name> <begin> <end>)')
    name, begin, end = fields
    if begin == '0':
        sys.exit('Error: please start counting with 1 in the blocks file')
    blocks.append([name, int(begin), int(end)])
    for line in blocksfile:
        fields = line.rstrip().split()
        if int(fields[1]) <= int(end):
            sys.exit('Error: overlapping blocks detected or blocks not correctly ordered')
        elif int(fields[1]) > int(end) + 1:
            i += 1
            blocks.append([f'unnamed{i}', int(end) + 1, int(fields[1]) - 1])
        name, begin, end = fields
        blocks.append([name, int(begin), int(end)])
    # Check if last one not more than length and if there is a chunk missing
    if int(end) > length:
        sys.exit('Error: blocks exceed the length of the alignment')
    elif int(end) < length:
        i += 1
        blocks.append([f'unnamed{i}', int(end) + 1, length])
    return blocks

# Options: step size of X% (default 1), single column, directly from original alignment
parser = argparse.ArgumentParser(description = "This script takes a multiple sequence alignment and removes compositionally biased sites.")
parser.add_argument("msa", help = "multiple sequence alignment")
parser.add_argument("output", help = "output file")
parser.add_argument("percentage", help = "remove X%% of alignment", type = int)
parser.add_argument("-f", metavar = "format", help = "alignment format (DEFAULT: FASTA)", default = "fasta")
parser.add_argument("-s", metavar = "stepsize", help = "remove sites by increments of X%% (DEFAULT: 1%%)", default = 1, type = int)
parser.add_argument("-b", metavar = "blocks", help = "file defining sections of the alignment (e.g., different markers) in format <name> <begin> <end>, fraction removed is reported for different blocks")
parser.add_argument("-i", help = "write all intermediate trimmed alignments for every step", action = "store_true")
parser.add_argument("-p", metavar = "prefix", help = "prefix of blocks and intermediate output files")
#parser.add_argument("-j", help = "create Jalview sequence features file", action = "store_true")
args = parser.parse_args()

# Parse alignment
align = AlignIO.read(args.msa, args.f)

# Define blocks in alignment
if args.b:
    with open(args.b) as blocksfile:
        blocks = parse_blocks(align, blocksfile)
    removed_per_block = {block[0]: [] for block in blocks}

# Calculate chi2 scores on full alignment
chi2, alphabet = calculate_chi2_scores(align)

# Calculate deltachi2 scores when removing each position and remove positions whose removal results in a lower chi2 score 
total_to_remove = args.percentage/100*align.get_alignment_length()
steps = args.percentage / args.s
increments = [int(i) for i in np.linspace(0, total_to_remove, int(steps)+1)]
trimmed_aln = align
original_positions = {i: i for i in range(1, align.get_alignment_length() + 1)}
for j, i in enumerate(range(len(increments) - 1)):
    no_positions = increments[i+1] - increments[i]
    if j > 0:
        chi2, _ = calculate_chi2_scores(trimmed_aln, report = False)
    print(f'Removing {no_positions} positions ({(j+1)*args.s}%)...')
    #position_deltachi2 = calculate_deltachi2_positions_old(trimmed_aln, chi2, alphabet)
    taxon_seqs = []
    total = ''
    for record in trimmed_aln:
        seq = str(record.seq)
        taxon_seqs.append(seq)
        total += seq.replace('-', '').replace('X', '').replace('*', '')
    freqs = {letter: total.count(letter) for letter in alphabet}
    position_deltachi2 = calculate_deltachi2_positions(trimmed_aln, freqs, len(total), taxon_seqs, chi2, alphabet)
    positions_to_remove = sorted(position_deltachi2, key = position_deltachi2.get)[:no_positions]
    positions_to_remove.sort()
    aln_pruned = trimmed_aln[:,:positions_to_remove[0]]
    previous_position = positions_to_remove[0]
    for position in positions_to_remove[1:]:
        aln_pruned += trimmed_aln[:,previous_position+1:position]
        previous_position = position
    aln_pruned += trimmed_aln[:,previous_position+1:]
    trimmed_aln = aln_pruned
    for original, position in original_positions.items():
        if position - 1 in positions_to_remove:
            original_positions[original] = 0
        else:
            change = sum([1 for pos in positions_to_remove if pos < position])
            original_positions[original] = position - change
    if args.b:
        for name, begin, end in blocks:
            removed = 0
            for position in range(begin, end+1):
                if original_positions[position] == 0:
                    removed += 1
            removed_per_block[name].append(removed)
    if i < len(increments) - 2 and args.i: # Not last iteration and intermediate alignments option
        out_name = f'removed_{j+1}pct.aln'
        if args.p:
            out_name = args.p + '.' + out_name
        AlignIO.write(trimmed_aln, out_name, args.f)

# Calculate chi2 scores on final trimmed alignment
print('\n\nCalculating chi2 scores on final alignment:')
calculate_chi2_scores(trimmed_aln, report = True)

# Write final trimmed alignment
AlignIO.write(trimmed_aln, args.output, args.f)

# Write positions removed per block
if args.b:
    if args.p:
        count_out_name = f'{args.p}.removed_per_block.tsv'
        fraction_out_name = f'{args.p}.removed_fraction_per_block.tsv'
    else:
        count_out_name = 'removed_per_block.tsv'
        fraction_out_name = 'removed_fraction_per_block.tsv'
    with open(count_out_name, 'w') as count_outfile, open(fraction_out_name, 'w') as fraction_outfile:
        steps = '\t'.join([f'{i}%' for i in range(args.s, args.percentage + 1, args.s)])
        count_outfile.write(f'Name\t{steps}\n')
        fraction_outfile.write(f'Name\t{steps}\n')
        for name, begin, end in blocks:
            block_size = end - begin + 1
            removed_counts = removed_per_block[name]
            out_string = '\t'.join([str(count) for count in removed_counts])
            count_outfile.write(f'{name}\t{out_string}\n')
            removed_fractions = '\t'.join([str(round(count / block_size, 4)) for count in removed_counts])
            fraction_outfile.write(f'{name}\t{removed_fractions}\n')