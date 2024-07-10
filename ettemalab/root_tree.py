#!/usr/bin/env python

import sys
import argparse
from ete3 import PhyloTree

# Parse arguments
parser = argparse.ArgumentParser(description = "This script roots a phylogenetic tree on the specified taxa.")
parser.add_argument("tree", help = "tree file")
parser.add_argument("outgroup", help = "comma-separated list of outgroup taxa")
parser.add_argument("-o", metavar = "output", help = "name of rooted tree file (DEFAULT: input + '.rooted')")
args = parser.parse_args()

tree = PhyloTree(args.tree, format = 1)
taxa = args.outgroup.split(',')
if len(taxa) < 2:
    sys.exit("Provide at least two outgroup taxon labels.")
outgroup = tree.get_common_ancestor(taxa)
if outgroup.is_root():
    sys.exit('Common ancestor of specified outgroup taxa is already the root. If unrooted, try providing complement taxa.')
tree.set_outgroup(outgroup)
if args.o:
    outfile = args.o
else:
    outfile = args.tree + '.rooted'
tree.write(format = 1, outfile = outfile)