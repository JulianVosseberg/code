#!/usr/bin/python3

# Functions for using the Eukarya database

import sys

def get_supergroups():
    supergroups_file = open("/home/julian/julian2/snel-clan-genomes/eukarya/5_supergroups.tsv")
    supergroups2 = {}
    supergroups5 = {}
    for line in supergroups_file:
        line = line.rstrip()
        line = line.split("\t")
        abbr = line[0]
        supergroup = line[1]
        supergroups5[abbr] = supergroup
        if supergroup == 'Obazoa' or supergroup == 'Amoebozoa':
            supergroups2[abbr] = 'Opimoda'
        else:
            supergroups2[abbr] = 'Diphoda'  
    supergroups_file.close()
    return supergroups2, supergroups5

def seq_info(species):
    seq_species_info = {}
    for spec in species:
        try:
            info_file = open('/hosts/linuxhome/memory/julian2/snel-clan-genomes/eukarya/proteomes_metadata/' + spec + '.metadata.txt')
        except IOError:
            print('Metadata file for %s could not be opened' % spec, file = sys.stderr)
            return False
        seq_species_info[spec] = {}
        for line in info_file:
            line = line.rstrip()
            fields = line.split('\t')
            if fields[7] != '1': # Not longest transcript
                continue
            seqid = fields[0]
            gene_symbol = fields[2]
            description = fields[4]
            seq_species_info[spec][seqid] = (gene_symbol, description)
        info_file.close()
    return seq_species_info
