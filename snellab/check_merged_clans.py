#!/usr/bin/python

## Script checks if Pfams that belong to distinct clans are merged

import sys; import re

clustering_file = sys.argv[1]
clans_file = sys.argv[2]
output_file = sys.argv[3]

try:
    clans_file = open(clans_file)
except IOError:
    print "%s cannot be opened" % clans_file
    sys.exit()

try:
    clustering_file = open(clustering_file)
except IOError:
    print "%s cannot be opened" % clustering_file
    sys.exit()

clans = {}
for line in clans_file:
    line = line.split("\t")
    pfam = line[0]
    clan = line[1]
    clans[pfam] = clan
clans_file.close()

output_file = open(output_file, "w")
output_file.write("Cluster\tNo_clans\tClans\n")

for line in clustering_file:
    pfams = re.findall("PF\d{5}", line)
    line = line.split("\t")
    cluster = line[0]
    merged_clans = []
    if len(pfams) > 1:
        for pfam in pfams:
            clan = clans[pfam]
            if clan != "" and clan not in merged_clans:
                merged_clans.append(clan)
        count = len(merged_clans)
        if count > 1:
            merged_string = "\t".join(merged_clans)
            output_file.write("\t".join([cluster, str(count), merged_string]) + "\n")

clustering_file.close()
output_file.close()
