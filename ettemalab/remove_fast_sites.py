#!/usr/bin/env python3

import argparse
from Bio import AlignIO
import pandas as pd
import sys

# Set arguments
parser = argparse.ArgumentParser(description = "This script removes (fast-evolving) sites from a multiple sequence alignment based on an inferred rates file. Note: resulting alignment as stdout.")
parser.add_argument("msa", help = "multiple sequence alignment (fasta)")
parser.add_argument("rates", help = "rates file inferred by IQ-TREE")
parser.add_argument("percentage", help = "remove X%% of sites", type = int)
parser.add_argument("-f", metavar = "format", help = "alignment format (DEFAULT: FASTA)", default = "fasta")
parser.add_argument("-s", help = "remove slowest instead of fastest sites", action = "store_true")
args = parser.parse_args()

# Parse alignment
align = AlignIO.read(args.msa, args.f)

# Parse rates and get indices of sites to remove
rate_df = pd.read_table(args.rates, skiprows = 8)
if args.s:
    sites_sorted = rate_df.Rate.sort_values(ascending = True).index
else:
    sites_sorted = rate_df.Rate.sort_values(ascending = False).index
number = round(len(sites_sorted) * args.percentage / 100)
to_remove = sorted(sites_sorted[:number])

# Remove sites
position = to_remove[0]
trimmed_aln = align[:, :position]
for new_position in to_remove[1:]:
    trimmed_aln += align[:, position+1:new_position]
    position = new_position

# Write output alignment
AlignIO.write(trimmed_aln, sys.stdout, args.f)