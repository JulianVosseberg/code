#!/usr/bin/python
import os; import sys; import re

## Script checks if the cluster file contains more than a specified threshold of sequences from at least two eukaryotic supergroups
## Otherwise, the cluster is moved to a separate directory

tax_file = sys.argv[1]
clusters_dir = sys.argv[2]
one_group_dir = sys.argv[3]
threshold = sys.argv[4]

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
unikonta = ["Amoebozoa", "Unknown", "Opisthokonta"]
sar = ["SAR", "Haptophyceae"]
archaepl = ["Archeaplastida", "Cryptophyta"]

for line in tax_file: # Specific for the eukarya v4 species file
    line = line.rstrip()
    line = line.split("\t")
    if line[0] == "core set":
        continue
    else:
        supergroup = line[3]
        if supergroup in unikonta:
            species_group[line[2]] = "Unikonta"
        elif supergroup in sar:
            species_group[line[2]] = "SAR"
        elif supergroup in archaepl:
            species_group[line[2]] = "Archaeplastida"
        elif supergroup == "Excavata":
            species_group[line[2]] = supergroup
        else:
            print "%s not recognised" % supergroup
tax_file.close()

for filename in os.listdir(clusters_dir):
    file_pattern = re.compile("[0-9]+\.fa$")
    if file_pattern.search(filename):
        file_path = clusters_dir + "/" + filename
	try:
            cluster_file = open(file_path)
        except IOError:
            print "%s cannot be opened" % filename
            sys.exit()
        supergroup_counter = {}
        for line in cluster_file:
            if line.startswith(">"):
                euk_pattern = re.compile("^>[A-Z]{4}")
                euk_match = re.search(euk_pattern, line)
                if euk_match:
                    euk_abbr = euk_match.group(0)[1:]
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
