#!/usr/bin/python
import sys; import re

## Script calculates for the PANTHER HMM based clustering the number of PANTHERs and COGs/ENOGs in a cluster

clustering_file = sys.argv[1]
output_file = sys.argv[2]

try:
    clustering_file = open(clustering_file)
except IOError:
    print "%s cannot be opened" % clustering_file
    sys.exit()

try:
    output_file = open(output_file, "w")
except IOError:
    print "%s cannot be opened" % output_file
    sys.exit()

output_file.write("Cluster\tNumber of PANTHERs\tNumber of COGs/ENOGs\tTotal\n")
    
for line in clustering_file:
    line = line.rstrip()
    cluster = line[0:line.find("\t")]
    panthers = line.count("\tPTHR")
    cogs_enogs = line.count("\tCOG") + line.count("\tENOG")
    total = line.count("\t")
    if total != panthers + cogs_enogs:
        print "For %s total does not correspond to the number of PANTHERs and COGs/ENOGs" % cluster
    output_file.write("\t".join([cluster, str(panthers), str(cogs_enogs), str(total)]) + "\n")

clustering_file.close()
output_file.close()
