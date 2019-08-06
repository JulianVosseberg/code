#!/usr/bin/env python3

import sys

pfams = sys.argv[1:]
out_files = {}
for pfam in pfams:
    out_files[pfam] = open(pfam + '.hhm', 'w')
pfam_hhm_file = open('/home/julian/julian2/pfam_hmm/hhm/pfam_hhm.ffdata')
found = False
curr_pfam = ''
for line in pfam_hhm_file:
    if not found:
        if line.startswith('NAME'):
            pfam = line[line.find('P'):line.find('.')]
            if pfam not in pfams:
                continue
            else:
                curr_pfam = pfam
                pfams.remove(pfam)
                found = True
                print('HHsearch 1.5', file = out_files[curr_pfam])
                print(line.rstrip(), file = out_files[curr_pfam])
    else:
        print(line.rstrip(), file = out_files[curr_pfam])
        if line.startswith('//'):
            found = False
            if len(pfams) == 0:
                break
else:
    print(', '.join(pfams) + ' not found', file = sys.stderr)
pfam_hhm_file.close()
for out_file in out_files.values():
    out_file.close()
