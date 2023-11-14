#!/usr/bin/env python

import sys
import glob
import os

infaa_dir = sys.argv[1]
blastfile = sys.argv[2]
outdir = sys.argv[3]

if not os.path.isdir(outdir):
    os.makedirs(outdir)

# Parse fasta files, create dictionaries with speciesIDs and sequenceIDs, write fasta files
speciesIDs = {}
sequenceIDs = {}
for species_count, faafilename in enumerate(glob.glob(infaa_dir + '/*.fa*')):
    speciesIDs[species_count] = os.path.basename(faafilename)
    with open(faafilename) as faafile, open(f'{outdir}/Species{species_count}.fa', 'w') as newfaa:
        seq_count = 0
        for line in faafile:
            if line[0] == '>':
                sequenceID = f'{species_count}_{seq_count}'
                newfaa.write(f'>{sequenceID}\n')
                sequenceIDs[sequenceID] = line.rstrip()[1:]
                seq_count += 1
            else:
                newfaa.write(line)
with open(outdir + '/SpeciesIDs.txt', 'w') as species_file:
    for speciesID, species in speciesIDs.items():
        species_file.write(f'{speciesID}: {species}\n')
original_newIDs = {}
with open(outdir + '/SequenceIDs.txt', 'w') as seq_file:
    for sequenceID, originalID in sequenceIDs.items():
        seq_file.write(f'{sequenceID}: {originalID}\n')
        if ' ' in originalID:
            originalID = originalID[:originalID.find(' ')]
        original_newIDs[originalID] = sequenceID

# Parse BLAST file
blast_results = {}
with open(blastfile) as blast_original:
    for line in blast_original:
        fields = line.split('\t')
        query = original_newIDs.get(fields[0], 'NA')
        subject = original_newIDs.get(fields[1], 'NA')
        if query == 'NA' or subject == 'NA':
            continue
        new_line = query + '\t' + subject + '\t' + '\t'.join(fields[2:])
        blast_search = query[:query.find('_')] + '_' + subject[:subject.find('_')]
        if blast_search not in blast_results:
            blast_results[blast_search] = []
        blast_results[blast_search].append(new_line)
for search, hits in blast_results.items():
    with open(f'{outdir}/Blast{search}.txt', 'w') as blast_out:
        blast_out.write(''.join(hits))
            
# SpeciesIDs.txt
# SequenceIDs.txt
# Species0.fa, Species1.fa, ...
# Blast0_0.txt, Blast0_1.txt, ...