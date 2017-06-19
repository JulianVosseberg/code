#!/usr/bin/python

import sys

def usage():
    print """\tUsage: parse_hmmsearch_Pfam.py <hmmsearch.domtblout> <output.tsv>\
\n\tThis script assigns sequences to their best Pfam hits."""

if len(sys.argv) != 3:
    usage(); sys.exit()

domtblout = sys.argv[1]
output_file = sys.argv[2]
overlap_cutoff = 15

try:
    hmmfile = open(domtblout)
except IOError:
    print domtblout, "could not be opened"; sys.exit()

seq_pfam = {}

for line in hmmfile:
    if line.startswith("#"):
        continue
    line = line.rstrip()
    line = line.split()
    seq = line[0]
    pfam = line[4]
    pfam = pfam[:pfam.find(".")]
    i_value = float(line[12])
    seq_start = int(line[19])
    seq_stop = int(line[20])
    if seq in seq_pfam: # Check for overlap, if there is: best hit, else: split sequence
        for index, hit in enumerate(seq_pfam[seq]):
            start_pr = hit[2]
            stop_pr = hit[3]
            overlap = max(seq_stop, stop_pr) - min(seq_start, start_pr) - ((stop_pr - start_pr) + (seq_stop - seq_start))
            if overlap <= -1 * overlap_cutoff: # Overlapping
                if i_value < hit[1]: # Check which hit is better
                    seq_pfam[seq][index] = [pfam, i_value, seq_start, seq_stop]
                break
        else: # Not overlapping with any other Pfams
            print seq, "hit by different non-overlapping Pfams and therefore assigned to multiple Pfams"
            seq_pfam[seq].append([pfam, i_value, seq_start, seq_stop])
    else:
        seq_pfam[seq] = [[pfam, i_value, seq_start, seq_stop]]
hmmfile.close()

pfam_seqs = {}
for seq in seq_pfam:
    for hit in seq_pfam[seq]:
        pfam = hit[0]
        start = hit[2]
        stop = hit[3]
        hit_info = seq + " (" + str(start) + ".." + str(stop) + ")"
        if pfam in pfam_seqs:
            pfam_seqs[pfam] += [hit_info]
        else:
            pfam_seqs[pfam] = [hit_info]

output_file = open(output_file, "w")         
for pfam in pfam_seqs:
    output_file.write(pfam + "\t" + "\t".join(pfam_seqs[pfam]) + "\n")
output_file.close()
