#!/usr/bin/python

import sys; import re

def usage():
    print "\tUsage: check_merged_clans_PTHR.py <clustering file> <COG Pfam hits> <Eukarya4 PANTHER hits> <Eukarya4 Pfam hits> [<clans file> <output_file>]"
    print "\n\tScript checks if in a cluster there are proteins hit by Pfams that belong to distinct clans"

no_args = len(sys.argv)
if no_args < 5 or no_args > 7:
    usage(); sys.exit()

clustering_file = sys.argv[1] #BBH_BH_clusters.PANTHER11.1-eggnog4.proteins.prokaryotes.domtblout_parsed_PANTHER.eukaryotic_update
cog_pfam_hits = sys.argv[2] #eggnog-prok_pfam.domtblout_parsed_EggNOG_no_coord
euk_pthr_hits = sys.argv[3] #eukarya_panther.domtblout.parsed.eukaryotic_memberships
euk_pfam_hits = sys.argv[4] #no_clust_parsed_Pfam
clans_file = "/hosts/linuxhome/memory/julian2/pfam_hmm/Pfam-A.clans.tsv"
output_file = "merged_clans.tsv"
if no_args == 6 or no_args == 7:
    clans_file = sys.argv[5]
if no_args == 7:
    output_file = sys.argv[6]

try:
    clans_file = open(clans_file)
except IOError:
    print clans_file, "cannot be opened"
    sys.exit()

clans = {}
for line in clans_file:
    line = line.split("\t")
    pfam = line[0]
    clan = line[1]
    clans[pfam] = clan
clans_file.close()

try:
    euk_pfam_hits = open(euk_pfam_hits)
except IOError:
    print euk_pfam_hits, "cannot be opened"
    sys.exit()

euk_pfam = {}
for line in euk_pfam_hits:
    line = line.rstrip()
    line = line.split("\t")
    pfam = line[0]
    for euk_seq in line[1:]:
        euk_pfam[euk_seq] = pfam
euk_pfam_hits.close()

try:
    euk_pthr_hits = open(euk_pthr_hits)
except IOError:
    print euk_pthr_hits, "cannot be opened"
    sys.exit()

pthr_pfams = {}
for line in euk_pthr_hits:
    line = line.rstrip()
    line = line.split("\t")
    pthr = line[0]
    pfams = []
    for euk_seq in line[1:]:
        if euk_seq in euk_pfam:
            pfam = euk_pfam[euk_seq]
            if pfam not in pfams:
                pfams.append(pfam)
        else:
            print euk_seq, "does not have a Pfam hit"
    pthr_pfams[pthr] = pfams
euk_pthr_hits.close()

try:
    cog_pfam_hits = open(cog_pfam_hits)
except IOError:
    print cog_pfam_hits, "cannot be opened"
    sys.exit()

cog_pfams = {}
for line in cog_pfam_hits:
    line = line.rstrip()
    line = line.split("\t")
    cogs = line[0].split("_")
    pfams = line[1:]
    for cog in cogs:
        cog_pfams[cog] = pfams
cog_pfam_hits.close()

try:
    clustering_file = open(clustering_file)
except IOError:
    print clustering_file, "cannot be opened"
    sys.exit()

output_file = open(output_file, "w")
output_file.write("Cluster\tNo_clans\tClans\n")

pthr_pattern = re.compile("PTHR")
cog_pattern = re.compile("COG")
enog_pattern = re.compile("ENOG")

for line in clustering_file:
    line = line.rstrip()
    line = line.split("\t")
    cluster = line[0]
    pfams_cluster = []
    for cog_pthr in line[1:]:
        if re.match(pthr_pattern, cog_pthr):
            if "_" in cog_pthr:
                print "In %s a split PANTHER (%s) found, number of clans found may be too high" % (cluster, cog_pthr)
                cog_pthr = re.sub("_\d+", "", cog_pthr)
            if cog_pthr in pthr_pfams:
                pfams = pthr_pfams[cog_pthr]
            else:
                print cog_pthr, "does not have any Pfam hits"
        elif re.match(cog_pattern, cog_pthr) or re.match(enog_pattern, cog_pthr):
            if "_" in cog_pthr:
                cogs = cog_pthr.split("_")
                pfams = []
                for cog in cogs:
                    if cog in cog_pfams:
                        pfams.extend(cog_pfams[cog])
                    else:
                        print cog, "does not have any Pfam hits"
            elif cog_pthr in cog_pfams:
                pfams = cog_pfams[cog_pthr]
            else:
                print cog_pthr, "does not have any Pfam hits"
        else:
            print cog_pthr, "not recognised"
        pfams_cluster.extend(pfams)
    pfams_cluster = set(pfams_cluster)
    clans_cluster = []
    for pfam in pfams_cluster:
        if "_" in pfam:
            print pfam, "was split during clustering"
            pfam = re.sub("_\d+", "", pfam)
        clan = clans[pfam]
        if clan not in clans_cluster:
            clans_cluster.append(clan)
    if len(clans_cluster) > 1:
        output_file.write("\t".join([cluster, str(len(clans_cluster))] + clans_cluster) + "\n")

output_file.close()
clustering_file.close()
