#!/usr/bin/env python

import numpy as np
import sys
import os
import argparse
from ete3 import PhyloTree

# Functions
def parse_metadata(metadatafile, cont_weight, complete):
    genome_info = {}
    with open(metadatafile) as metadata:
        for line in metadata:
            fields = line.rstrip().split('\t')
            name = fields[0]
            compl = float(fields[2])
            cont = float(fields[3])
            if complete and fields[4] == 'yes':
                gqs = 101.0
            else:
                gqs = compl - cont_weight * cont
            genome_info[name] = {'markers': int(fields[1]), 'completeness': compl, 'contamination': cont, 'GQS': gqs}
    return genome_info

def annotate_node(node, genome_info):
    if node.name in genome_info:
        node.add_features(gqs=genome_info[node.name]['GQS'], completeness=genome_info[node.name]['completeness'], contamination=genome_info[node.name]['contamination'], markers=genome_info[node.name]['markers'])

def calculate_RED_values(tree, genome_info):
    tree.add_feature("red", 0)
    for node in tree.iter_descendants():
        if node.is_leaf():
            red = 1
            annotate_node(node, genome_info)
        else:
            parent = node.up
            p = parent.red
            d = node.get_distance(parent)
            u = np.mean([parent.get_distance(leaf) for leaf in node])
            red = p + d/u*(1-p)
            # RED as internal node name to display it in figtree
            node.name = red
        node.add_feature("red", red)

# Parse arguments
parser = argparse.ArgumentParser(description = "This script takes a rooted phylogenetic tree, calculates RED values and selects the best representatives.")
parser.add_argument("tree", help = "rooted tree file")
parser.add_argument("metadata", help = "tab-delimited file with leaf names, number of markers, completeness estimate, contamination estimate and if it is a closed genome (yes/no)")
parser.add_argument("-o", metavar = "outdir", help = "directory for output files (DEFAULT: current)")
parser.add_argument("-r", metavar = "RED_value", help = "RED value used to define clades to pick representatives from (DEFAULT: 0.8)", default = 0.8, type = float)
parser.add_argument("-i", help = "only calculate RED values, no selection (DEFAULT: off)", action = "store_true")
parser.add_argument("-c", help = "preferentially select a complete genome (DEFAULT: off)", action = "store_true")
parser.add_argument("-w", metavar = "weight", help = "weight given to the contamination estimate (genome quality score = completeness - X * contamination)", default = 5, type = int)
group = parser.add_mutually_exclusive_group()
group.add_argument("-a", metavar = "leaf1,leaf2,leaf3", help = "only select representatives from the subtree that contains the defined leaves (DEFAULT: off)")
group.add_argument("-s", metavar = "subgroup", help = "only select representatives from a clade whose leaf labels start with this prefix (DEFAULT: off)")
parser.add_argument("-m", metavar = "markers", help = "only select genomes with at least X markers (DEFAULT: 30)", default = 30, type = int)
parser.add_argument("-t", metavar = "contamination", help = "only select genomes with at most X%% estimated contamination (DEFAULT: 10.0)", default = 10.0, type = float)
parser.add_argument("-p", metavar = "completeness", help = "only select genomes with at least X%% estimated completeness (DEFAULT: 50.0)", default = 50.0, type = float)
args = parser.parse_args()

# Parse metadata and tree
if args.o:
    outdir = args.o
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
else:
    outdir = '.'
genome_info = parse_metadata(args.metadata, args.w, args.c)
tree = PhyloTree(args.tree)
calculate_RED_values(tree, genome_info)

# # If only calculating RED values, write tree with RED values and stop here
if args.i:
    tree.write(outfile = f"{outdir}/tree_with_RED_values.nw", format = 1)
    sys.exit()

# Subtree
if args.a:
    tree = tree.get_common_ancestor(args.a.split(','))
elif args.s:
    tree = tree.get_common_ancestor([leaf.name for leaf in tree if leaf.name.startswith(args.s)])

# Pick representatives
selection = []
tree.add_feature('collapse', 'no')
for node in tree.iter_descendants():
    if node.up.collapse in ('yes', 'here'):
        node.add_feature('collapse', 'yes')
    elif node.red > args.r:
        node.add_feature('collapse', 'here')
        best_genome = ''
        best_gqs = 0
        for leaf in node.iter_leaves():
            if leaf.markers < args.m:
                continue
            if leaf.contamination >= args.t or leaf.completeness < args.p:
                continue
            if leaf.gqs > best_gqs:
                best_genome = leaf.name
                best_gqs = leaf.gqs
        if best_genome == '':
            sys.stderr.write(f'Warning: no genome fulfills requirements in clade {",".join([leaf.name for leaf in node.iter_leaves()])}.\n')
        else:
            selection.append(best_genome)
    else:
        node.add_feature('collapse', 'no')

# Output: list with representatives and pruned tree (only subgroup)
with open(f'{outdir}/selected_representatives.list', 'w') as repr_list:
    repr_list.write('\n'.join(selection) + '\n')
tree.prune(selection, preserve_branch_length = True)
tree.write(outfile = f"{outdir}/selected_representatives.nw")
