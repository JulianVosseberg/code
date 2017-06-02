#!/usr/bin/python
import os; import sys; import re

## Script checks if there are eukaryotic sequences in a cluster, if not it is moved to a prok_only directory
## Assumption: Fasta headers of eukaryotic sequences start with abbreviation (as is the case for Eukarya v4, whereas prokaryotic sequences start with the taxonomy ID)

clusters_dir = sys.argv[1]
prok_only_dir = sys.argv[2]

try:
    os.makedirs(prok_only_dir)
except OSError:
    print "%s cannot be created" % prok_only_dir
    sys.exit()

for filename in os.listdir(clusters_dir):
    file_pattern = re.compile("^C[0-9]+\.fa")
    if file_pattern.match(filename):
        try:
            cluster_file = open(filename)
        except IOError:
            print "%s cannot be opened" % filename
            sys.exit()
        euk_found = False
        for line in cluster_file:
            if line.startswith('>'):
                euk_pattern = re.compile("^>[A-Z]{4}")
                if euk_pattern.match(line):
                    euk_found = True
                    break
        cluster_file.close()
        if not euk_found:
            os.rename(filename, (prok_only_dir + "/" + filename))
