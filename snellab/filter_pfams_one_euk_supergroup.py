#!/usr/bin/python
import os; import sys; import re

## Script checks if the cluster file contains more than a specified threshold of sequences from at least two eukaryotic supergroups
## Otherwise, the cluster is moved to a separate directory

tax_file = '/home/julian/julian2/snel-clan-genomes/eukarya/5_supergroups.tsv'
pfams_dir = sys.argv[1]
one_group_dir = sys.argv[2]
threshold = sys.argv[3]

try:
    os.makedirs(one_group_dir)
except OSError:
    print "%s cannot be created" % one_group_dir
    sys.exit()

try:
    tax_file = open(tax_file)
except IOError:
    print "%s cannot be opened" % tax_file
    sys.exit()

species_group = {}
opimoda = ["Obazoa", "Amoebozoa"]

for line in tax_file: # Specific for the eukarya v4 species file
    line = line.rstrip()
    line = line.split("\t")
    supergroup = line[1]
    if supergroup in opimoda:
        species_group[line[0]] = "Opimoda"
    else:
        species_group[line[0]] = supergroup
tax_file.close()

for filename in os.listdir(pfams_dir):
    file_pattern = re.compile("\.fa$")
    if file_pattern.search(filename):
        file_path = pfams_dir + "/" + filename
	try:
            cluster_file = open(file_path)
        except IOError:
            print "%s cannot be opened" % filename
            sys.exit()
        supergroup_counter = {}
        for line in cluster_file:
            if line.startswith(">"):
                euk_abbr = line[1:5]
                if euk_abbr not in species_group:
                    print "%s not found" % euk_abbr
                else:
                    supergroup = species_group[euk_abbr]
                    if supergroup in supergroup_counter:
                        supergroup_counter[supergroup] += 1
                    else:
                        supergroup_counter[supergroup] = 1
        cluster_file.close()
        no_groups = 0
        for supergroup in supergroup_counter:
            if supergroup_counter[supergroup] >= int(threshold):
                no_groups += 1
        if no_groups < 2:
            os.rename(file_path, (one_group_dir + "/" + filename))
