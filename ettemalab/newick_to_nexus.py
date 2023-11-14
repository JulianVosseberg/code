#!/usr/bin/env python

from Bio import Phylo
import sys
import os

tree_file = sys.argv[1]
Phylo.convert(tree_file,  'newick', tree_file + '.tmp.nex', 'nexus')
with open(tree_file + '.tmp.nex') as tmp_nexus, open(tree_file + '.nex', 'w') as nexus:
    for line in tmp_nexus:
        if 'TaxLabels' in line:
            line = line.replace(';', '')
            nexus.write(line.replace(' ', '\n ') + ' ;\n')
        else:
            nexus.write(line)
os.remove(tree_file + '.tmp.nex')