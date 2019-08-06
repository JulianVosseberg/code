#!/usr/bin/python

import sys

# New version of parse_hmm_search_Pfam.py. This one does not extract the best hits based on the order of hits.

def usage():
    print """\tUsage: parse_hmmsearch_Pfam_2.py <hmmsearch.domtblout> <output.tsv>\
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

seq_hits = {}

for line in hmmfile:
    if line.startswith("#"):
        continue
    line = line.rstrip()
    line = line.split()
    seq = line[0]
    pfam = line[4]
    pfam = pfam[:pfam.find(".")]
    i_value = float(line[12])
    hmm_start = int(line[15])
    hmm_stop = int(line[16])
    seq_start = int(line[19])
    seq_stop = int(line[20])
    try:
        seq_hits[seq].append([pfam, i_value, seq_start, seq_stop, hmm_start, hmm_stop])
    except KeyError:
        seq_hits[seq] = [[pfam, i_value, seq_start, seq_stop, hmm_start, hmm_stop]]
hmmfile.close()

pfam_seqs = {}
for seq, hits in seq_hits.items():
    if len(hits) > 1: # Check for remaining overlapping Pfams
        indices_to_remove = []
        for x in xrange(len(hits) - 1):
            for y in xrange(1, len(hits)):
                if x == y:
                    continue
                overlap = max(hits[x][3], hits[y][3]) - min(hits[x][2], hits[y][2]) - ((hits[x][3] - hits[x][2]) + (hits[y][3] - hits[y][2]))
                if overlap <= -1 * overlap_cutoff: # Overlap
                    if hits[x][1] < hits[y][1]: # x is better
                        indices_to_remove.append(y)
                    else:
                        indices_to_remove.append(x)
        if len(indices_to_remove) != 0:
            for index in set(indices_to_remove):
                print 'overlap', seq, hits[index][0]
            indices_to_keep = set(range(len(hits))) - set(indices_to_remove)
            hits = [hits[index] for index in indices_to_keep]
        hits.sort(key = lambda hit: hit[2]) # Sort hits based on hit start coordinates
    pfams = [hit[0] for hit in hits]
    if len(pfams) != len(set(pfams)): # So, duplicates
        pfams_seen = []
        pfams_dupl = []
        for pfam in pfams:
            if pfam in pfams_seen:
                if pfam not in pfams_dupl:
                    pfams_dupl.append(pfam)
            else:
                pfams_seen.append(pfam)
        for pfam in pfams_dupl:
            pfam_hits = [hit[2:6] + [index] for index, hit in enumerate(hits) if hit[0] == pfam]
            for i in xrange(len(pfam_hits) - 1):
                overlap = max(pfam_hits[i][3], pfam_hits[i+1][3]) - min(pfam_hits[i][2], pfam_hits[i+1][2]) - ((pfam_hits[i][3] - pfam_hits[i][2]) + (pfam_hits[i+1][3] - pfam_hits[i+1][2]))
                if overlap > -1 * overlap_cutoff: # Non-overlapping HMMs
                    #print 'split Pfam' # Check if orientation is okay (only then merge, otherwise history not clear, e.g. Ras)
                    if pfam_hits[i][2] < pfam_hits[i+1][2]:
                        print 'insert', seq, hits[pfam_hits[i][-1]][0]
                        hits[pfam_hits[i][-1]].append('c')
    for hit in hits:
        pfam = hit[0]
        start = hit[2]
        stop = hit[3]
        if hit[-1] == 'c':
            hit_info = seq + " c(" + str(start) + ".." + str(stop) + ")" # put 'c' before brackets to indicate merge with next hits wished
        else:
            hit_info = seq + " (" + str(start) + ".." + str(stop) + ")"
        try:
            pfam_seqs[pfam].append(hit_info)
        except KeyError:
            pfam_seqs[pfam] = [hit_info]

output_file = open(output_file, "w")
for pfam in pfam_seqs:
    output_file.write(pfam + "\t" + "\t".join(pfam_seqs[pfam]) + "\n")
output_file.close()
