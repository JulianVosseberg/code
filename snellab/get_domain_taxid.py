#!/usr/bin/python
import sys; import re; import os

## Script counts the number of bacterial, archaeal and eukaryotic sequences in each cluster file, based on the taxonomy ID of the prokaryotic sequences (eukaryotic ones are characterised by an abbreviation i/o tax ID)

clusters_dir = sys.argv[1]
taxid_file = sys.argv[2]
output_file = sys.argv[3]

try:
    taxid_file = open(taxid_file)
except IOError:
    print "%s cannot be opened" % taxid_file
    sys.exit()

prok_domains = {}
for line in taxid_file:
    line = line.rstrip()
    line = line.split("\t")
    domain = line[0]
    if domain == "A" or domain == "B":
        prok_domains[line[2]] = domain
taxid_file.close()

output_file = open(output_file, "w")
output_file.write("Cluster\tBacteria\tArchaea\tEukaryotes\tNot recognised\n")

for filename in os.listdir(clusters_dir):
    file_pattern = re.compile("^C[0-9]+\.fa")
    if file_pattern.match(filename):
        try:
            cluster_file = open(clusters_dir + "/" + filename)
        except IOError:
            print "%s cannot be opened" % filename
            sys.exit()
        domain_counter = {"A":0, "B":0, "E":0, "N":0}
        for line in cluster_file:
            if line.startswith('>'):
                euk_pattern = re.compile("^>[A-Z]{4}")
                euk_match = re.search(euk_pattern, line)
                if euk_match:
                    domain_counter["E"] += 1
                else:
                    prok_pattern = re.compile("^>\d+\.")
                    prok_match = re.search(prok_pattern, line)
                    if prok_match:
                        tax_id = prok_match.group(0)[1:-1]
                        if tax_id not in prok_domains.keys():
                            print "Taxonomy ID %s in %s not recognised" % (tax_id, filename)
                            domain_counter["N"] += 1
                        else:
                            domain = prok_domains[tax_id]
                            domain_counter[domain] += 1
                    else:
                        print "%s in %s does not correspond to a bacterial/archaeal/eukaryotic sequence" % (line, filename)
        cluster_file.close()
        result = [filename[:-3], str(domain_counter["B"]), str(domain_counter["A"]), str(domain_counter["E"]), str(domain_counter["N"])]
        output_file.write("\t".join(result) + "\n")

output_file.close()
        
