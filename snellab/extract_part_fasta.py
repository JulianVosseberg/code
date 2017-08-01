#!/usr/bin/python

import sys

def usage():
    print "Usage: extract_part_fasta.py <fasta> <start> <end>"
    print "Script extracts for each sequence the part from the given start up to and including the given end cordinates and provides it as standard output"

if len(sys.argv) != 4:
    usage(); sys.exit()

try:
    fasta_file = open(sys.argv[1])
except IOError:
    print sys.argv[1], "cannot be opened"; sys.exit()

start = int(sys.argv[2]) - 1
end = int(sys.argv[3])

seq = ""
for line in fasta_file:
    line = line.rstrip()
    if line[0] == ">":
        print seq[start:end]
        print line
        seq = ""
    else:
        seq += line
else:
    print seq[start:end]
fasta_file.close()
