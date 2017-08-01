#!/usr/bin/python

import sys

def usage():
    print "Usage: number_multiple_domains.py <file.fa>"
    sys.exit()

if len(sys.argv) != 2:
    usage(); sys.exit()

fasta=sys.argv[1]

try:
    fastafile = open(fasta)
except IOError:
    print fasta, "cannot be opened"; sys.exit()

# Loop over file to get sequence IDs
seqids = []
multiple_counts = {}
for line in fastafile:
    if line.startswith(">"):
        line = line.rstrip()
        seqid = line[1:]
        if seqid not in seqids:
            seqids.append(seqid)
        else:
            if seqid not in multiple_counts:
                multiple_counts[seqid] = 0

# Back to start of file
fastafile.seek(0)

# Open output file
outputfile = open(fasta[:-3] + "_domains_numbered.fa", "w")
for line in fastafile:
    if line.startswith(">"):
        seqid = line.rstrip()[1:]
        if seqid in multiple_counts:
            number = multiple_counts[seqid] + 1
            outputfile.write(">" + seqid + "_domain_" + str(number) + "\n")
            multiple_counts[seqid] += 1
        else:
            outputfile.write(line)
    else:
        outputfile.write(line)

# Close files
outputfile.close()
fastafile.close()
