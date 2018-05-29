#! /usr/bin/env python3

import sys

def usage():
    sys.exit("\tUsage: get_alignment_length.py <fasta>\n\
    \tThis script calculates the length of an alignment in fasta format by calculating the length of the first sequence.")

if len(sys.argv) != 2:
    usage()

fasta_name = sys.argv[1]

try:
    fasta = open(fasta_name)
except IOError:
    sys.exit(fasta_name, 'could not be opened')

first = True
length = 0

for line in fasta:
    if first:
        if line.startswith('>'):
            first = False
        else:
            fasta.close()
            print('Error: file not in fasta format', file = sys.stderr); usage()
    else:
        if line.startswith('>'):
            print(fasta_name, length)
            fasta.close()
            break
        else:
            line = line.rstrip()
            length += len(line)
else:
    print('Only one sequence in %s with length %d' % (fasta_name, length))
