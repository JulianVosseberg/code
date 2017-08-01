#!/usr/bin/python

import sys

def usage():
	print "Usage: count_supergroups.py <supergroups.tsv> <species_abbreviation1> <species_abbreviation2>..."
	print "Provide a supergroups file in the format ABBREVIATION\\tSUPERGROUP and a list with the species abbreviations"

if len(sys.argv) < 3:
	usage(); sys.exit()

try:
	supergroups_file = open(sys.argv[1])
except IOError:
	print sys.argv[1], "cannot be opened"; sys.exit()

supergroups = {}
species_list = sys.argv[2:]
supergroups_counter = {}

for line in supergroups_file:
	line = line.rstrip()
	line = line.split("\t")
	abbr = line[0]
	supergroup = line[1]
	supergroups[abbr] = supergroup

supergroups_file.close()

counter = 0
for species in species_list:
	if species not in supergroups:
		print species, "could not be found in the supergroups file and not taken into account"
	else:
		counter += 1
		supergroup = supergroups[species]
		if supergroup not in supergroups_counter:
			supergroups_counter[supergroup] = 1
		else:
			supergroups_counter[supergroup] += 1

for supergroup in supergroups_counter:
	print supergroup, ":", supergroups_counter[supergroup]
print "Total:", counter
