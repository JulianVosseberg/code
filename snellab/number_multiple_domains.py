#!/usr/bin/env python2

import sys
import getopt
import re

def usage():
    sys.exit("Usage: number_multiple_domains.py [ -d <delimiter> ] <file.fa>\n\nIf a delimiter is given, the part behind this in the header is ignored")

optlist, args = getopt.getopt(sys.argv[1:], "d:")
delim = ''
if len(optlist) == 1:
    delim = optlist[0][1]
if len(args) != 1:
   usage()

fasta=args[0]

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
        if delim != '':
            m = re.search(delim, seqid)
            if m:
                seqid = line[1 : m.start() + 1]
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
        line = line.rstrip()
        seqid = line[1:]
        seqid_full = seqid
        if delim != '':
            m = re.search(delim, seqid)
            if m:
                seqid = line[1 : m.start() + 1]
                seqid_full = line[1:]
        if seqid in multiple_counts:
            number = multiple_counts[seqid] + 1
            outputfile.write(">" + seqid_full + "_domain_" + str(number) + "\n")
            multiple_counts[seqid] += 1
        else:
            outputfile.write(line + '\n')
    else:
        outputfile.write(line)

# Close files
outputfile.close()
fastafile.close()
