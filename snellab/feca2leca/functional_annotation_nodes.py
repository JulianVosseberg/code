#!/usr/bin/env python3

import sys
import re

def get_annotation(seqs, seqs_function, group):
    if group not in ('Function', 'Localisation'):
        sys.stderr.write('Error: annotation group not recognised.\n')
        return
    unkn_return = {'Function' : 'S', 'Localisation' : 'NA'}
    function_count = {}
    for seq in seqs:
        if '_' in seq:
            seq = seq[:seq.find('_')]
        try:
            functions = seqs_function[seq][group]
        except KeyError:
            sys.stderr.write(f'Warning: {seq} has no annotation for {group}.\n')
            continue
        for function in functions:
            function_count[function] = function_count.get(function, 0) + 1
    if len(function_count) > 1:
        sorted_function_count = sorted(function_count, key = function_count.get, reverse = True)
        if function_count[sorted_function_count[0]] > function_count[sorted_function_count[1]]:
            return sorted_function_count[0]
        else:
            return unkn_return[group]
    elif len(function_count) == 0:
        return unkn_return[group]
    else: # Just one function
        return function

interest_ids = {
'GO:0005576' : 'extracellular region',
'GO:0005618' : 'cell wall',
'GO:0005886' : 'plasma membrane',
'GO:0005730' : 'nucleolus',
'GO:0005654' : 'nucleoplasm',
'GO:0005635' : 'nuclear envelope',
'GO:0031410' : 'cytoplasmic vesicle',
'GO:0000228' : 'nuclear chromosome',
'GO:0005783' : 'endoplasmic reticulum',
'GO:0005794' : 'Golgi apparatus',
'GO:0005773' : 'vacuole',
'GO:0005777' : 'peroxisome',
'GO:0005929' : 'cilium',
'GO:0005768' : 'endosome',
'GO:0005829' : 'cytosol',
'GO:0005739' : 'mitochondrion',
'GO:0005856' : 'cytoskeleton'
}

# Load eggnogmapper output
seqs_function = {}
with open('/home/julian/julian2/eggnog4/eukarya_4.0.2/eukarya_4.0.2.emapper.annotations') as eggnog_annotation_file:
    for line in eggnog_annotation_file:
        line = line.rstrip()
        fields = line.split('\t')
        if len(fields) == 13:
            seq = fields[0]
            functions = fields[11].split(', ')
            seqs_function[seq] = {'Function': functions}
            go_terms = fields[5]
            seqs_function[seq]['Localisation'] = [go_term for go_term in go_terms.split(',') if go_term in interest_ids]

og_function = {}
with open('/home/julian/julian2/timing_dupl/4_phylogeny_pfam_pipeline2/1_bbhs_2_groups/4_ete_analysis/d20_l15/lecas.tsv') as lecas_file:
    lecas_file.readline()
    for line in lecas_file:
        line = line.rstrip()
        fields = line.split('\t')
        pfam = fields[0]
        og = fields[3]
        if pfam not in og_function:
            og_function[pfam] = {}
        seqs = fields[9].split(',')
        og_function[pfam][og] = {'Function' : get_annotation(seqs, seqs_function, 'Function')}
        og_function[pfam][og]['Localisation'] = get_annotation(seqs, seqs_function, 'Localisation')

# Integrate the localisation annotation for the LECAs with the duplication information
prok_euk_pattern = re.compile('(PF\d{5}_FECA[\d_]+)_(PF.+)')
with open('/home/julian/julian2/timing_dupl/4_phylogeny_pfam_pipeline2/1_bbhs_2_groups/5_assign_euk_clan/d20_l15/merged_pfams.tsv') as duplication_info_file:
    duplication_info_file.readline()
    for line in duplication_info_file:
        line = line.rstrip()
        fields = line.split('\t')
        name = fields[0]
        ancestry = fields[1]
        if '(maj)' in ancestry:
            ancestry = ancestry[:-6]
        euks = []
        duplicated = int(fields[2]) > 1
        if 'FECA' in name:
            # Eukaryote-only assigned to donation
            if '_PF' in name:
                match = re.search(prok_euk_pattern, name)
                if not match:
                    sys.stderr.write(f'Error in regex search for {name}!\n')
                    continue
                euks = match.group(2).split('_')
                feca = match.group(1)
            # Only donation, might still contain multiple FECAs (merged FECA)
            else:
                feca = name
            pfam = name[:7]
            feca_nos = feca[12:].split('_')
            for og in og_function[pfam]:
                if og[2:og.find('.')] in feca_nos:
                    og_function[pfam][og]['Ancestry'] = ancestry
                    og_function[pfam][og]['Duplicated'] = duplicated
        # Invention
        else:
            euks = name.split('_')
        for pfam in euks:
            for og in og_function[pfam]:
                og_function[pfam][og]['Ancestry'] = ancestry
                og_function[pfam][og]['Duplicated'] = duplicated

with open('/home/julian/julian2/timing_dupl/4_phylogeny_pfam_pipeline2/1_bbhs_2_groups/6_functional_annotation/lecas.tsv', 'w') as lecas_file:
    lecas_file.write('Pfam\tOG\tAncestry\tDuplicated\tFunction\tLocalisation\n')
    for pfam, ogs in og_function.items():
        for og, og_info in ogs.items():
            lecas_file.write(f'{pfam}\t{og}\t{og_info["Ancestry"]}\t{og_info["Duplicated"]}\t{og_info["Function"]}\t{interest_ids.get(og_info["Localisation"], "NA")}\n')

# Function duplication nodes
with open(f'/home/julian/julian2/timing_dupl/4_phylogeny_pfam_pipeline2/1_bbhs_2_groups/4_ete_analysis/d20_l15/duplication_lengths.tsv') as dupl_file:
    dupl_annotation = {}
    dupl_file.readline()
    for line in dupl_file:
        fields = line.split('\t')
        pfam = fields[0]
        if pfam not in dupl_annotation:
            dupl_annotation[pfam] = {}
        dupl = fields[3]
        daughters = fields[7].split(' - ')
        dupl_annotation[pfam][dupl] = {'Daughters' : daughters}
        try:
            dupl_annotation[pfam][dupl]['Dl'] = float(fields[9])
        except ValueError: # For PF00005 branch lengths removed
            dupl_annotation[pfam][dupl]['Dl'] = 'NA'
        functions1 = set([og_function[pfam][og]['Function'] for og in daughters[0].split(',')])
        functions2 = set([og_function[pfam][og]['Function'] for og in daughters[1].split(',')])
        ancestral = functions1.intersection(functions2)
        if ancestral == {'S'}:
            dupl_annotation[pfam][dupl]['Function'] = 'S'
        elif len(ancestral - {'S'}) == 1:
            #print(f"Ancestral function for {dupl}: {tuple(ancestral)[0]}")
            dupl_annotation[pfam][dupl]['Function'] = tuple(ancestral - {'S'})[0]
        elif len(ancestral) > 1:
            #sys.stderr.write(f'Warning: multiple ancestral functions for {dupl}.\n')
            dupl_annotation[pfam][dupl]['Function'] = 'S'
        else:
            # Here comes the tricky part!
            # If parent duplication and the function of this duplication node overlaps with one of the daugher functions...
            if len(dupl_annotation[pfam]) > 1: # There are parent nodes (not necessarily this FECA)
                # Combine all daughters of node of interest
                daughters_combined = set(daughters[0].split(',') + daughters[1].split(','))
                for parent, info in dupl_annotation[pfam].items():
                    if parent == dupl:
                        continue
                    parent_daughters = info['Daughters']
                    # If one of the parent daughters fully contains the daughters of the current node, it is the parent node
                    if set(parent_daughters[0].split(',')) == daughters_combined or set(parent_daughters[1].split(',')) == daughters_combined:
                        parent_function = info['Function']
                        if parent_function == 'S':
                            #sys.stderr.write(f'Warning: function of parent node of {dupl} unclear.\n')
                            dupl_annotation[pfam][dupl]['Function'] = 'S'
                        elif parent_function in functions1.union(functions2):
                            #print(f'Ancestral function for {dupl} (based on parent node): {parent_function}')
                            dupl_annotation[pfam][dupl]['Function'] = parent_function
                        else:
                            dupl_annotation[pfam][dupl]['Function'] = 'S'
                        break
                else:
                    #sys.stderr.write(f'Warning: no parent nodes for {dupl}, all other FECAs.\n')
                    dupl_annotation[pfam][dupl]['Function'] = 'S'
            else:
                #sys.stderr.write(f'Warning: ancestral function for {dupl} unclear.\n')
                dupl_annotation[pfam][dupl]['Function'] = 'S'
        go1 = set([og_function[pfam][og]['Localisation'] for og in daughters[0].split(',')])
        go2 = set([og_function[pfam][og]['Localisation'] for og in daughters[1].split(',')])
        ancestral = go1.intersection(go2)
        if ancestral == {'NA'}:
            dupl_annotation[pfam][dupl]['Localisation'] = 'NA'
        elif len(ancestral - {'NA'}) == 1:
            #print(f"Ancestral localisation for {dupl}: {slim_ids[tuple(ancestral)[0]]}")
            dupl_annotation[pfam][dupl]['Localisation'] = tuple(ancestral)[0]
        elif len(ancestral) > 1:
            #sys.stderr.write(f'Warning: multiple ancestral localisations for {dupl}.\n')
            dupl_annotation[pfam][dupl]['Localisation'] = 'NA'
        else:
            # Here comes the tricky part!
            # If parent duplication and the function of this duplication node overlaps with one of the daugher functions...
            if len(dupl_annotation[pfam]) > 1: # There are parent nodes (not necessarily this FECA)
                # Combine all daughters of node of interest
                daughters_combined = set(daughters[0].split(',') + daughters[1].split(','))
                for parent, info in dupl_annotation[pfam].items():
                    if parent == dupl:
                        continue
                    parent_daughters = info['Daughters']
                    # If one of the parent daughters fully contains the daughters of the current node, it is the parent node
                    if set(parent_daughters[0].split(',')) == daughters_combined or set(parent_daughters[1].split(',')) == daughters_combined:
                        parent_go = info['Localisation']
                        if parent_go == 'NA':
                            #sys.stderr.write(f'Warning: localisation of parent node of {dupl} unclear.\n')
                            dupl_annotation[pfam][dupl]['Localisation'] = 'NA'
                        elif parent_go in go1.union(go2):
                            #print(f'Ancestral localisation for {dupl} (based on parent node): {slim_ids[parent_go]}')
                            dupl_annotation[pfam][dupl]['Localisation'] = parent_go
                        else:
                            dupl_annotation[pfam][dupl]['Localisation'] = 'NA'
                        break
                else:
                    #sys.stderr.write(f'Warning: no parent nodes for {dupl}, all other FECAs.\n')
                    dupl_annotation[pfam][dupl]['Localisation'] = 'NA'
            else:
                #sys.stderr.write(f'Warning: ancestral localisation for {dupl} unclear.\n')
                dupl_annotation[pfam][dupl]['Localisation'] = 'NA'

with open('/home/julian/julian2/timing_dupl/4_phylogeny_pfam_pipeline2/1_bbhs_2_groups/6_functional_annotation/duplications.tsv', 'w') as dupl_output:
    dupl_output.write('Pfam\tDuplication\tFunction\tLocalisation\tDl\n')
    for pfam, duplications in dupl_annotation.items():
        for duplication, info in duplications.items():
            print(pfam, duplication, info['Function'], interest_ids.get(info['Localisation'], info['Localisation']), info['Dl'], sep = '\t', file = dupl_output)
