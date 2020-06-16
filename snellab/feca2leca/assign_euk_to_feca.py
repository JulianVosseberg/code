#!/usr/bin/env python3

# First part
import os
from subprocess import PIPE, Popen, call
import sys
import re

groups = sys.argv[1]
if groups == '2':
    path = '/home/julian/julian2/timing_dupl/4_phylogeny_pfam_pipeline2/1_bbhs_2_groups/4_ete_analysis/d20_l15/'
elif groups == '5':
    path = '/home/julian/julian2/timing_dupl/4_phylogeny_pfam_pipeline2/2_bbhs_5_groups/4_ete_analysis/d20_l15/'
elif groups == 'r2':
    path = '/home/julian/julian2/timing_dupl/6_nee_revision/1_kclust_incl_asgard/5_ete_analysis/bbhs_2_groups/'
else:
    sys.exit("Please provide a valid argument: 2 (2 groups) or 5 (5 groups)")

# Characterise clans and do simple merge if possible
clans_file = open('/home/julian/julian2/pfam_hmm/Pfam-A.clans.tsv')
clan_pfams = {}
pfam_clan = {}
for line in clans_file:
    line = line.rstrip()
    fields = line.split('\t')
    pfam = fields[0]
    clan = fields[1]
    if clan == '':
        continue
    clan_name = fields[2]
    pfam_clan[pfam] = clan
    clan_pfams[clan] = [clan_name, [], []]
clans_file.close()

pe_pfams = []
e_pfams = []
pfam_fecas = {}
pfam_lecas = {}
pfam_unknowns = {}
leca_pfams = open(path + 'stats_per_pfam.tsv')
for line in leca_pfams:
    if line.startswith('Pfam'):
        continue
    line = line.rstrip()
    fields = line.split('\t')
    lecas = int(fields[4])
    if lecas == 0: # At least 1 LECA
        continue
    pfam = fields[0]
    if pfam not in pfam_clan:
        continue
    pfam_unknowns[pfam] = int(fields[5])
    if fields[1] == 'NA': # FECAs for euk are NA
        e_pfams.append(pfam)
        pfam_lecas[pfam] = lecas
    else:
        pe_pfams.append(pfam)
        mfecas = int(fields[3])
        pfam_fecas[pfam] = mfecas
leca_pfams.close()

for pfam, clan in pfam_clan.items():
    if pfam in pe_pfams:
        index = 1
    elif pfam in e_pfams:
        index = 2
    else:
        continue
    clan_pfams[clan][index].append(pfam)

if groups == '2':
    assign_path = '/home/julian/julian2/timing_dupl/4_phylogeny_pfam_pipeline2/1_bbhs_2_groups/5_assign_euk_clan/d20_l15/'
elif groups == 'r2':
    assign_path = '/home/julian/julian2/timing_dupl/6_nee_revision/1_kclust_incl_asgard/6_assign_euk_clan/bbhs_2_groups/'
else:
    assign_path = '/home/julian/julian2/timing_dupl/4_phylogeny_pfam_pipeline2/2_bbhs_5_groups/5_assign_euk_clan/d20_l15/'
outfile = open(assign_path + 'clan_pfams.tsv', 'w')
print('Clan\tName\tProk+euk\tEuk', file = outfile)
mergedfile = open(assign_path + 'merged_pfams.tsv', 'w')
print('Merged name\tAncestry\tLECAs\tUnknowns\tMethod', file = mergedfile)
pfams_merged = []
fecas_merged = []
to_merge = {}
for clan, content in clan_pfams.items():
    if len(content[2]) == 0: # No euk Pfams
        continue
    elif len(content[2]) == 1:
        if len(content[1]) == 0: # Only one euk Pfam
             continue
    print(clan, content[0], ','.join(content[1]), ','.join(content[2]), sep = '\t', file = outfile)
    if len(content[1]) == 0:
        pfams = content[2]
        pfams_merged.extend(pfams)
        name = '_'.join(pfams)
        lecas = sum([pfam_lecas[pfam] for pfam in pfams])
        unknowns = sum([pfam_unknowns[pfam] for pfam in pfams])
        print(name, 'Eukaryotic', lecas, unknowns, 'simple', sep = '\t', file = mergedfile)
    elif len(content[1]) == 1:
        pe_pfam = content[1][0]
        if pfam_fecas[pe_pfam] == 1:
            info = Popen(f'grep -m 1 {pe_pfam} {path}/fecas_strict.tsv', stdout = PIPE, shell = True).communicate()[0]
            info = info.decode().split('\t')
            fecaname = info[1]
            e_pfams = content[2]
            name = '_'.join([pe_pfam, fecaname] + e_pfams)
            pfams_merged.extend(e_pfams)
            fecas_merged.append(pe_pfam + '_' + fecaname)
            lecas = int(info[4]) + sum([pfam_lecas[pfam] for pfam in e_pfams])
            print(name, info[3], lecas, 'TBD', 'simple', sep = '\t', file = mergedfile)
        else:
            to_merge[clan] = {'pe' : content[1], 'e' : content[2]}
    else:
        to_merge[clan] = {'pe' : content[1], 'e' : content[2]}
outfile.close()

merged_feca = {}
for clan, pfams in to_merge.items():
    os.makedirs(clan)
    sys.stderr.write(f'Selecting sequences per FECA for {clan}...\n')
    for pfam in pfams['pe']:
        call(f'get_feca_sequence_ids.py {pfam} {groups} {clan}', shell = True)
    sys.stderr.write('Finished selecting sequences!\n')
    sys.stderr.write(f'Performing hhsearch for {clan}...\n')
    call(f'assign_euk_to_feca.sh {groups} {clan} {",".join(pfams["e"])}', shell = True)
    sys.stderr.write(f'Finished hhsearch for {clan}!\n')
    for e_pfam in pfams['e']:
        best_hit = Popen(f'grep -m 1 " 1 PF" {assign_path}{clan}/euk_only/{e_pfam}.hhr', stdout = PIPE, shell = True).communicate()[0]
        feca = best_hit.decode().split()[1] #re.search('PF\d{5}_FECA\w+', best_hit).group(0)
        try:
            merged_feca[feca].append(e_pfam)
        except KeyError:
            merged_feca[feca] = [e_pfam]
for feca, e_pfams in merged_feca.items():
    fecas_merged.append(feca)
    pfams_merged.extend(e_pfams)
    name = '_'.join([feca] + e_pfams)
    feca_search = feca.replace('_FECA', '\tFECA')
    info = Popen(f'grep -m 1 "{feca_search}" {path}/fecas_strict.tsv', stdout = PIPE, shell = True).communicate()[0]
    info = info.decode().split('\t')
    lecas = int(info[4]) + sum([pfam_lecas[pfam] for pfam in e_pfams])
    print(name, info[3], lecas, 'TBD', 'search', sep = '\t', file = mergedfile)

# Combine with other FECAs and inventions
statsfile = open(path + 'stats_per_pfam.tsv')
for line in statsfile:
    if line.startswith('Pfam'):
        continue
    if 'NA' in line:
        fields = line.split('\t')
        e_pfam = fields[0]
        if e_pfam in pfams_merged:
            continue
        print(e_pfam, 'Eukaryotic', fields[4], fields[5], 'NA', sep = '\t', file = mergedfile)
fecas_file = open(path + 'fecas_strict.tsv')
for line in fecas_file:
    if line.startswith('Pfam'):
        continue
    fields = line.split('\t')
    feca = fields[0] + '_' + fields[1]
    if feca in fecas_merged:
        continue
    print(feca, fields[3], fields[4], 'TBD', 'NA', sep = '\t', file = mergedfile)
mergedfile.close()
