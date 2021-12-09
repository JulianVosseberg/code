#!/usr/bin/env python2

import sys; import os

def usage():
    print """\tUsage: merge_sequences_per_Pfam.py <pfams.list> <directory 1> (<directory 2> ...)\
\n\tThis script reads the list of hit Pfams and merges the corresponding sequence files."""

if len(sys.argv) < 3:
    usage(); sys.exit()

pfams = []
try:
    pfams_list = open(sys.argv[1])
except IOError:
    print sys.argv[1], "not found"; sys.exit()
for pfam in pfams_list:
    pfam = pfam.rstrip()
    pfams.append(pfam)
pfams_list.close()

try:
    os.mkdir("combined")
except OSError:
    print "Directory 'combined' cannot be created"; sys.exit()

for pfam in pfams:
    filenames = []
    for directory in sys.argv[2:]:
        if not os.path.isdir(directory):
            print directory, "not found"; sys.exit()
        filename = directory + "/" + pfam + ".fa"
        if os.path.isfile(filename):
            filenames.append(filename)
    if len(filenames) == 0:
        print pfam, "not found in any provided directory"
    else:
        os.system("cat " + " ".join(filenames) + " > combined/" + pfam + ".fa")
