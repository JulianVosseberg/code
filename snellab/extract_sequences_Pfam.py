#!/usr/bin/python

import sys; import getopt

def usage():
    print """\tUsage: extract_sequences_Pfam.py [-d] <tsv file 1> <fasta sequences 1> (<tsv file 2> ...)\
\n\tThis script reads the sequence IDs for each hit Pfam and extracts the corresponding sequence entirely (default) or only the hit domain (-d argument)."""

optlist, args = getopt.getopt(sys.argv[1:], "d")
domains = False
if len(optlist) == 1:
    domains = True
if len(args) % 2 != 0 or len(args) == 0:
    usage(); sys.exit()

pfams = []

for j,i in enumerate(xrange(0, len(args), 2)):
    try:
        tsv = open(args[i])
    except IOError:
        print args[i], "not found"; sys.exit()
    try:
        fasta = open(args[i+1])
    except IOError:
        print args[i+1], "not found"; sys.exit()
    seqs = {}; seq_id = ""
    for line in fasta:
        line = line.rstrip()
        if line.startswith(">"):
            seq_id = line[1:]
            seqs[seq_id] = ""
        else:
            seqs[seq_id] += line
    fasta.close()
    for line in tsv:
        line = line.rstrip()
        line = line.split("\t")
        pfam = line[0]
        output_file = open(pfam + "_" + str(j) + ".fa", "w")
        if pfam not in pfams:
            pfams.append(pfam)
        for seq in line[1:]:
            seq_info = seq.split()
            seq_id = seq_info[0]
            if seq_id not in seqs:
                print seq, "not found. Script aborted."; sys.exit()
            if domains:
                coordinates = seq_info[1].split("..")
                start = int(coordinates[0][1:])
                stop = int(coordinates[1][:-1])
                output_file.write(">" + seq_id + "\n" + seqs[seq_id][start-1:stop] + "\n")
            else:
                output_file.write(">" + seq_id + "\n" + seqs[seq_id] + "\n")
        output_file.close()
    tsv.close()

pfam_list = open("pfams.list", "w")
pfam_list.write("\n".join(pfams))
