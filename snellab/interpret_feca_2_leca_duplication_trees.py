#!/home/julian/.local/anaconda3/bin/python3

# Load modules
from ete3 import PhyloTree
from ete3 import NCBITaxa
from numpy import median
from eukarya import *
import re
import sys
import getopt
import os
import collections
ncbi = NCBITaxa()

# Usage
def usage():
    print('\n\tUsage: interpret_feca_2_leca_duplication_trees.py -t <tree.contree> [ -p <prefix> ] [ -o <output_dir> ] [ -b <bbhs.tsv> ] [ -e ] [ -f ] [ -i ] [ -r ] [ -h ]', file = sys.stderr)
    print('\nThis script identifies FECA-2-LECA duplications, determines the best prokaryotic outgroup (if any), and performs a branch length analysis', file = sys.stderr)
    print('\nOptions:\n\t-t: tree file\n\t-p: prefix for output files (DEFAULT: basename tree file)\n\t-o: directory for output files (DEFAULT: current)\
    \n\t-b: use list of BBHs for OG annotation\n\t-e: only eukaryotes (DEFAULT: off)\
    \n\t-f: use only farthest leaf for rooting (DEFAULT: off)\
    \n\t-i: filter interspersing prokaryotes (DEFAULT: off)\
    \n\t-r: use information of which other sequences are represented by a ScrollSaw sequence\
    \n\t-h: help', file = sys.stderr)
    sys.exit()

# Distinguish eukaryotic and prokaryotic leaves
def annotate_prokaryotic_eukaryotic_leaves(tree, euk_only):
    if euk_only:
        for leaf in tree:
            taxid = leaf.name[0:4]
            leaf.add_features(taxid = taxid, supergroup2 = supergroups2[taxid], supergroup5 = supergroups5[taxid])
    else:
        for leaf in tree:
            if re.match('^\d', leaf.name):
                taxid = int(leaf.name[:leaf.name.find('.')])
                lineage = ncbi.get_lineage(taxid)
                leaf.add_features(taxid = taxid, sp = ncbi.get_taxid_translator([taxid])[taxid], prok_euk = 'Prokaryote')
            else:
                taxid = leaf.name[0:4]
                leaf.add_features(taxid = taxid, supergroup2 = supergroups2[taxid], supergroup5 = supergroups5[taxid], prok_euk = 'Eukaryote')

# Get monophyletic eukaryotic clades
def get_euk_clades(tree):
    root = tree.search_nodes(prok_euk = 'Prokaryote')[0] # Root on first prokaryotic sequence
    tree.set_outgroup(root)
    clades = tree.get_monophyletic(values = ['Eukaryote'], target_attr = 'prok_euk')
    return clades

# Determine if clade fulfills FECA-2-LECA criteria: both Opimoda and Diphoda present
def feca2leca(leaves):
    opimoda = 0
    diphoda = 0
    for leaf in leaves:
        if leaf.supergroup2 == 'Opimoda':
            opimoda += 1
        elif leaf.supergroup2 == 'Diphoda':
            diphoda += 1
        if opimoda > 0 and diphoda > 0:
            break
    else:
        return False
    return True

# Determine if a node fulfills the FECA-2-LECA duplication criteria: both daughters have both Opimoda and Diphoda
def duplication_check(leaves1, leaves2):
    if feca2leca(leaves1) and feca2leca(leaves2):
        return True
    else:
        return False

# Reroot on farthest leaf
def reroot(euk_clade, tree):
    tree.set_outgroup(euk_clade) # Root on this eukaryotic clade
    sister = euk_clade.get_sisters()[0]
    farthest = sister.get_farthest_leaf()[0]
    tree.set_outgroup(farthest) # Root on the leaf farthest from this eukaryotic clade (can be a false positive for example)
    return farthest.name

# Determine both possible prokaryotic sister groups in an unrooted way or a rooted way using the rooting on the farthest leaf
def get_prokaryotic_sister(euk_clade, tree, farthest):
    if farthest:
        sister = euk_clade.get_sisters()[0] # Should be checked if there are any eukaryotic sequences in the sister group
        prok_leaves_sister = sister.search_nodes(prok_euk = 'Prokaryote')
        if len(tree) - (len(euk_clade) + len(sister)) == 1: # So only 1 other non-sister sequence
            support = 'NA'
        else:
            support = euk_clade.up.support
        prok_taxids = []
        for prok in prok_leaves_sister:
            prok_taxids.append(prok.taxid)
        sp_tree = ncbi.get_topology(prok_taxids) # Get NCBI species tree, to get the identity of the LCA
        lca = sp_tree.taxid
        lca_name = ncbi.translate_to_names([lca])[0]
        return lca, lca_name, support
    else:    
        sister = euk_clade.get_sisters()[0] # Should be checked if there are any eukaryotic sequences in the sister group
        prok_leaves_sister = sister.search_nodes(prok_euk = 'Prokaryote')
        other_prok_leaves = set(tree.search_nodes(prok_euk = 'Prokaryote')) - set(prok_leaves_sister)
        lcas = []
        if len(tree) - (len(euk_clade) + len(sister)) == 1: # So only 1 other non-sister sequence
            supports = ['NA']
        else:
            supports = [euk_clade.up.support]
        if len(sister) == 1:
                supports.append('NA')
        else:
            supports.append(sister.support)
        for i, group in enumerate([prok_leaves_sister, other_prok_leaves]):
            prok_taxids = []
            for prok in group: # Collect tax ids of prokaryotic sister leaves
                prok_taxids.append(prok.taxid)
            sp_tree = ncbi.get_topology(prok_taxids) # Get NCBI species tree, to get the identity of the LCA
            lca = sp_tree.taxid
            if lca == 1224: # Proteobacteria
                for proteo in sp_tree:
                    lineage = ncbi.get_lineage(proteo.name)
                    if not 28211 in lineage and not 1236 in lineage and not 28216 in lineage: # So, not an alpha/gamma/beta
                        lca_name = 'Proteobacteria'
                        break
                else: # So, only aerobic proteobacteria
                    lca = 'aerprot'
                    lca_name = 'Aerobic proteobacteria'
            elif lca == 2157: # Archaea
                for arch in sp_tree:
                    lineage = ncbi.get_lineage(arch.name)
                    if not 1935183 in lineage and not 1783275 in lineage: # So, not an Asgard or TACK
                        lca_name = 'Archaea'
                        break
                else: # So, only Asgards + TACK
                    lca = 'asgtack'
                    lca_name = 'Asgard+TACK supergroup'
            else:
                lca_name = ncbi.translate_to_names([lca])[0]
            support = supports[i]
            lcas.append((lca, lca_name, support))
        return lcas

# Classify the prokaryotic sister-group
def classify_sister(lca): # To do: add aerobic proteo superclass and TACK+Asgard supersuperphylum
    if lca == 'aerprot':
        return 'aerobic proteobacterial'
    elif lca == 'asgtack':
        return 'Asgardian+TACK'
    ancestors = ncbi.get_lineage(lca)
    if 28211 in ancestors:
        return 'alphaproteobacterial'
    elif 1935183 in ancestors:
        return 'Asgardian'
    elif 28216 in ancestors:
        return 'betaproteobacterial'
    elif 1236 in ancestors:
        return 'gammaproteobacterial'
    elif 1783275 in ancestors:
        return 'TACK'
    else:
        ranks = ncbi.get_rank(ancestors)
        names = ncbi.get_taxid_translator(ancestors)
        for taxon in ranks: # Return phylum and if that is not present, then lowest group
            if ranks[taxon] == 'phylum':
                return names[taxon]
        else:
            return names[ancestors[-1]]
            #print('Error: prokaryotic sister group', lca, 'not identified', file = sys.stderr)

# Annotate eukaryotic nodes in a euk-only tree
def annotate_duplications_midpoint_rooted(tree):
    for node in tree.traverse("preorder"):
        if node.is_leaf():
            continue
        parent = node.up
        leaves1 = node.children[0].get_leaves()
        leaves2 = node.children[1].get_leaves()
        if not node.is_root():
            if parent.identity == 'LECA' or parent.identity == 'post-LECA':
                node.add_features(identity = 'post-LECA')
                continue
        if duplication_check(leaves1, leaves2):
            node.add_features(identity = "duplication")
        elif feca2leca(node.get_leaves()):
            node.add_features(identity = 'LECA')
        else:
            node.add_features(identity = 'unknown')

# Annotate eukaryotic nodes
def annotate_eukaryotic_nodes(clades):
    for ca in clades:
        if not feca2leca(ca.get_leaves()):
            continue
        for node in ca.traverse('preorder'):
            if node.is_leaf():
                continue
            leaves1 = node.children[0].get_leaves()
            leaves2 = node.children[1].get_leaves()
            if duplication_check(leaves1, leaves2):
                node.add_features(identity = 'duplication')
            elif feca2leca(node.get_leaves()):
                node.add_features(identity = 'LECA')
            else:
                node.add_features(identity = 'post-LECA')
        for leca in ca.iter_search_nodes(identity = 'LECA'):
            if len(leca.search_nodes(identity = 'LECA')) > 1: # So, if there is another potential LECA further down, identity changed to unknown
                if len(leca.search_nodes(identity = 'duplication')) > 0: # So, duplication nodes downstream
                    leca.identity = 'unknown'
                else: # Only LECAs downstream
                    for other_leca in leca.search_nodes(identity = 'LECA')[1:]:
                        other_leca.identity = 'post-LECA(ch)'
#     Print per FECA, number of duplications, LECAs and unknowns to output file

# Get non-FECA sister
def get_non_feca_sister(non_feca_node, tree): # In a rooted way (on farthest leaf)
    tree.set_outgroup(non_feca_node)
    farthest_leaf = non_feca_node.get_sisters()[0].get_farthest_leaf()[0]
    tree.set_outgroup(farthest_leaf)
    sister = non_feca_node.get_sisters()[0]
    prok_leaves_sister = sister.search_nodes(prok_euk = 'Prokaryote')
    support = non_feca_node.up.support
    prok_taxids = []
    for prok in prok_leaves_sister:
        prok_taxids.append(prok.taxid)
    sp_tree = ncbi.get_topology(prok_taxids) # Get NCBI species tree, to get the identity of the LCA
    lca = sp_tree.taxid
    lca_name = ncbi.translate_to_names([lca])[0]
    fecas = []
    for feca_clade in sister.iter_search_nodes(feca = True):
        fecas.append(feca_clade.feca_no)
    if len(fecas) > 0:
        print('Warning: for non-FECA eukaryotic sequences in %s the sister group does contain a FECA' % prefix, file = sys.stderr)
        lca_name += '+' + '+'.join(fecas) 
    return lca, lca_name, support

# Calculate branch lengths for a single FECA
def calculate_branch_lengths(tree, lecas, lepca, sister):
    prok_branch_lengths = []
    for prok in sister:
        prok_branch_lengths.append(tree.get_distance(lepca, prok))
    prok_bl_med = median(prok_branch_lengths)
    raw_stem_lengths = []
    stem_lengths = []
    ebls = []
    for leca in lecas:
        euk_branch_lengths = []
        for euk in leca:
            euk_branch_lengths.append(tree.get_distance(leca, euk))
        euk_bl_med = median(euk_branch_lengths)
        rsl = tree.get_distance(lepca, leca)
        sl = rsl / euk_bl_med
        ebls.append(euk_bl_med)
        raw_stem_lengths.append(rsl)
        stem_lengths.append(sl)
    rsl_med = median(raw_stem_lengths)
    sl_med = median(stem_lengths)
    ebls_med = median(ebls)
    return prok_bl_med, rsl_med, sl_med, ebls_med

# Calculate duplication lengths --> inconsistent!!!
def calculate_median_duplication_length(tree, duplication, lecas):
    raw_dupl_lengths = []
    dupl_lengths = []
    for leca in lecas:
        euk_branch_lengths = []
        for euk in leca:
            euk_branch_lengths.append(tree.get_distance(leca, euk))
        ebl_med = median(euk_branch_lengths)
        raw_dupl_length = tree.get_distance(leca, duplication)
        raw_dupl_lengths.append(raw_dupl_length)
        dupl_length = raw_dupl_length / ebl_med
        dupl_lengths.append(dupl_length)
    rdls_med = median(raw_dupl_lengths)
    dls_med = median(dupl_lengths)
    return rdls_med, dls_med

def prok_filter(tree):
    # Note: if there is only one prokaryotic leaf, this one will be removed.
    euk_leave_names = [leaf.name for leaf in tree.iter_search_nodes(prok_euk = 'Eukaryote')]
    tree.set_outgroup(euk_leave_names[0])
    filtered_leaves = []
    for clade in tree.get_monophyletic(values = ['Prokaryote'], target_attr = 'prok_euk'):
        if len(clade) != 1: # So, only singletons
            continue
        parent = clade.up
        prok_name = clade.name
        sister = clade.get_sisters()[0]
        if parent.up.is_root(): # Immediate sister is where the tree is rooted on
            for child in sister.get_children():
                if len(child.search_nodes(prok_euk = 'Prokaryote')) == 0:
                    print('Interspersing prokaryote type 1:', clade.name, file = sys.stderr)
                    filtered_leaves.append(clade)
                    break
        elif len(sister.search_nodes(prok_euk = 'Prokaryote')) != 0:
            tree.set_outgroup(parent)
            clade = tree&prok_name
            parent = clade.up
            sister = clade.get_sisters()[0]
            new_sister = parent.get_sisters()[0]
            if len(new_sister.search_nodes(prok_euk = 'Prokaryote')) == 0:
                for child in sister.get_children():
                    if len(child.search_nodes(prok_euk = 'Prokaryote')) == 0:
                        print('Interspersing prokaryote type 2:', clade.name, file = sys.stderr)
                        filtered_leaves.append(clade)
                        tree.set_outgroup(tree&euk_leave_names[0])
                        break
                else:
                    tree.set_outgroup(tree&euk_leave_names[0])
            else:
                tree.set_outgroup(tree&euk_leave_names[0])
        else:
            p_sister = parent.get_sisters()[0]
            if len(p_sister.search_nodes(prok_euk = 'Prokaryote')) == 0:
                print('Interspersing prokaryote type 3:', clade.name, file = sys.stderr)
                filtered_leaves.append(clade)
            else:
                tree.set_outgroup(p_sister)
                clade = tree&prok_name
                parent = clade.up
                new_sister = parent.get_sisters()[0]
                if len(new_sister.search_nodes(prok_euk = 'Prokaryote')) == 0:
                    print('Interspersing prokaryote type 4:', clade.name, file = sys.stderr)
                    filtered_leaves.append(clade)
                tree.set_outgroup(tree&euk_leave_names[0])
    if len(filtered_leaves) == 0:
        pruned = False
        return(pruned)
    leaves_kept = [leaf for leaf in tree if leaf not in filtered_leaves]
    tree.prune(leaves_kept, preserve_branch_length = True)
    pruned = True
    return(pruned)

def reannotate_bbhs(clade, bbh_pairs):
    leca_nodes = clade.search_nodes(identity = 'LECA')
    unknown_nodes = clade.search_nodes(identity = 'unknown')
    pot_ogs = {}
    unknown_leaves = []
    for i,leca in enumerate(leca_nodes):
        leca.add_features(pot_og = 'pOG' + str(i))
        leave_names = [leaf.name for leaf in leca.iter_leaves()]
        pot_ogs[i] = leave_names
    for unknown in unknown_nodes:
        children = unknown.get_children()
        for child in children:
            if child.is_leaf() or child.identity == 'post-LECA':
                leave_names = [leaf.name for leaf in child.iter_leaves()]
                break
        else:
            print('Error in the annotation', file = og_annotation)
        unknown_leaves.extend(leave_names)
    for node in clade.traverse():
        node.del_feature('identity') # Erase previous identity
        node.add_features(leca_og = False)
    comb_ogs = []
    unknown_ogs = {}
    for i in pot_ogs:
        unknown_ogs[i] = []
    for pair in bbh_pairs:
        og_bbhs = []
        unknown = []
        for seq in pair:
            for pot_og in pot_ogs:
                if seq in pot_ogs[pot_og]:
                    og_bbhs.append(pot_og)
                    break
            else: # Either unknown or non-FECA
                unknown.append(seq)
        if len(og_bbhs) == 2:
            og1, og2 = og_bbhs
            if og1 != og2:
                print('Potential OGs %d and %d should be fused' % (og1, og2), file = og_annotation)
                new_og_set = set(og_bbhs)
                og_sets_to_replace = []
                for index,og_set in enumerate(comb_ogs):
                    if og1 in og_set or og2 in og_set: # Overlap; one can connect two previous sets, so all have to be checked
                        og_sets_to_replace.append(og_set)
                        new_og_set = og_set | new_og_set # Combine with rest
                for og_set in og_sets_to_replace:
                    comb_ogs.remove(og_set)
                comb_ogs.append(new_og_set)
        elif len(og_bbhs) == 1:
            unknown_seq = unknown[0]
            og = og_bbhs[0]
            if unknown_seq in unknown_leaves:
                print('Unclassified sequence %s should be assigned to OG %d' %(unknown_seq, og), file = og_annotation)
                unknown_ogs[og].append(unknown_seq)
            else:
                print('Warning: sequence %s not in this FECA clade should be assigned to OG %d' %(unknown_seq, og), file = og_annotation)
        else:
            seq1, seq2 = unknown
            if seq1 in unknown_leaves and seq2 in unknown_leaves:
                print('BBH between unclassified sequences %s and %s' %(seq1, seq2), file = og_annotation)
            elif seq1 in unknown_leaves or seq2 in unknown_leaves:
                if seq1 in unknown_leaves:
                    unknown_seq = seq1; other = seq2
                else:
                    unknown_seq = seq2; other = seq1
                    print('Unclassified sequences %s best hit by sequence %s not in this FECA clade' %(unknown_seq, other), file = og_annotation)
    pot2final_ogs = {}
    number_pogs = list(range(len(pot_ogs)))
    for pot_og in number_pogs:
        for og_set in comb_ogs:
            if pot_og in og_set:
                pot2final_ogs[pot_og] = list(og_set)
                for other_og in og_set - set([pot_og]):
                    number_pogs.remove(other_og)
                break
        else:
            pot2final_ogs[pot_og] = [pot_og]
    final_ogs = {}
    for og in pot2final_ogs:
        sequences = []
        for pot_og in pot2final_ogs[og]:
            sequences_pog = pot_ogs[pot_og] + unknown_ogs[pot_og]
            sequences.extend(sequences_pog)
        final_ogs[og] = sequences
    for og in final_ogs:
        og_node = clade.get_common_ancestor(final_ogs[og])
        og_node.add_features(leca_og = True, fin_og = 'fOG' + str(og))
    for og_node in clade.iter_search_nodes(leca_og = True):
        og_nodes_down = og_node.search_nodes(leca_og = True)
        if len(og_nodes_down) > 1:
            for og_node_down in og_nodes_down[1:]:
                print('%s not monophyletic and therefore %s merged into this fOG' %(og_node.fin_og, og_node_down.fin_og), file = og_annotation)
                og_node_down.add_features(leca_og = False, previous_fog = og_node_down.fin_og)
                og_node_down.del_feature('fin_og')
    for og_node in clade.iter_search_nodes(leca_og = True):
        moving_up = True
        while(moving_up):
            if og_node == clade: # Stop going up if FECA is reached
                break
            sister = og_node.get_sisters()[0]
            if len(sister.search_nodes(leca_og = True)) != 0: # Stop if the sister contains a LECA
                moving_up = False
            else:
                fin_og = og_node.fin_og
                og_node.add_features(leca_og = False)
                og_node.del_feature('fin_og')
                leaves_added = [leaf.name for leaf in sister.iter_leaves()]
                print('Included %s in %s' %(','.join(leaves_added), fin_og), file = og_annotation)
                og_node = og_node.up
                og_node.add_features(leca_og = True, fin_og = fin_og)
    for node in clade.traverse():
        if node.is_leaf():
            continue
        if node.leca_og:
            node.add_features(identity = 'LECA')
            continue
        children = node.get_children()
        if len(children[0].search_nodes(leca_og = True)) > 0 and len(children[1].search_nodes(leca_og = True)) > 0:
            node.add_features(identity = 'duplication')
        elif len(children[0].search_nodes(leca_og = True)) > 0 or len(children[1].search_nodes(leca_og = True)) > 0:
            node.add_features(identity = 'unknown')

def assign_all_seqs(tree, all_seqs, euk_only):
    if euk_only:
        tree_seqs = [leaf.name for leaf in tree]
    else:
        tree_seqs = [leaf.name for leaf in tree.iter_search_nodes(prok_euk = 'Eukaryote')]
    other_seqs = set(all_seqs) - set(tree_seqs)
    human_represent = {}
    for seq in tree_seqs:
        if 'HSAP' in seq:
            human_represent[seq] = seq # Human seq represented by itself
    if euk_only:
        blast = open('/hosts/linuxhome/memory/julian2/pfam_hmm/domains_extracted/combined/t1/euk_only/blast_output/' + prefix + '/' + prefix + '_blastp.txt')
    else:
        blast = open('/hosts/linuxhome/memory/julian2/pfam_hmm/domains_extracted/combined/t1/blast_output/' + prefix + '/' + prefix + '_blastp.txt')
    best_blast_hits = {}
    for line in blast:
        line = line.rstrip()
        fields = line.split('\t')
        query = fields[0]
        if query in tree_seqs:
            continue
        hit = fields[1]
        if hit not in tree_seqs:
            continue
        bit_score = float(fields[11])
        if query not in best_blast_hits:
            best_blast_hits[query] = (hit, bit_score)
        else:
            best_score = best_blast_hits[query][1]
            if bit_score > best_score:
                best_blast_hits[query] = (hit, bit_score)
    blast.close()
    representing = {}
    for query in best_blast_hits:
        best_hit = best_blast_hits[query][0]
        if best_hit in representing:
            representing[best_hit].append(query)
        else:
            representing[best_hit] = [query]
        if 'HSAP' in query:
            human_represent[query] = best_hit
    return representing, human_represent

def infer_coverage_redundancy(node, representing):
    number_species = {'Obazoa':123, 'Amoebozoa':6, 'RASH':43, 'Archaeplastida':24, 'Excavata':13}
    all_species = []
    species_counts = {}
    supergroup_counts = {'Obazoa':0, 'Amoebozoa':0, 'RASH':0, 'Archaeplastida':0, 'Excavata':0}
    group_species_counts = {'Obazoa':[], 'Amoebozoa':[], 'RASH':[], 'Archaeplastida':[], 'Excavata':[]}
    for leaf in node: # node can also be a list of leaves instead
        name = leaf.name
        sp = leaf.taxid
        species = [seq[0:4] for seq in representing.get(name,'')]
        species.append(sp)
        all_species.extend(species)
    for spec in all_species:
        if spec in species_counts:
            species_counts[spec] += 1
        else:
            supergroup = supergroups5[spec]
            supergroup_counts[supergroup] += 1
            species_counts[spec] = 1
    coverages = []
    for group, counts in supergroup_counts.items():
        coverages.append(counts / number_species[group])
    coverage_av = sum(coverages) / 5
    for spec, counts in species_counts.items():
        group = supergroups5[spec]
        group_species_counts[group].append(counts)
    redundancy_med = median([median(counts) for counts in group_species_counts.values() if len(counts) > 0])
    return coverage_av, redundancy_med, list(species_counts.keys())

def annotate_overlap_all_assigned(feca, overlap_criterion = 10, coverage_criterion = 0.15):
    for node in feca.traverse('preorder'):
        coverage, copies, repr_species = infer_coverage_redundancy(node, representing)
        node.add_features(coverage = coverage, redundancy = copies, repr_species = repr_species, identity = '?')
    for node in feca.traverse('preorder'):
        if node.is_leaf():
            continue
        leaves = node.get_leaves()
        if not feca2leca(leaves):
            continue
        daughters = node.get_children()
        species_overlap = set(daughters[0].repr_species) & set(daughters[1].repr_species)
        if len(species_overlap) >= overlap_criterion:
            if daughters[0].coverage >= coverage_criterion and daughters[1].coverage >= coverage_criterion:
                if duplication_check(daughters[0].get_leaves(), daughters[1].get_leaves()):
                    node.add_features(identity = 'duplication', overlap = len(species_overlap))
                    continue
        if node == feca:
            if node.coverage >= coverage_criterion:
                node.add_features(identity = 'LECA')
            else:
                node.add_features(coverage = coverage)
                return False
        elif node.up.identity == 'duplication':
            node.add_features(identity = 'LECA')
        else:
            node.add_features(identity = 'post-LECA')
    for leca in feca.iter_search_nodes(identity = 'LECA'):
        if len(leca.search_nodes(identity = 'duplication')) > 0:
            leca.identity = 'unknown'
        else:
            lecas_down = leca.search_nodes(identity = 'LECA')
            if len(lecas_down) > 1:
                for leca_down in lecas_down[1:]:
                    leca_down.identity = '?'
    return True

def duplication_check_unrooted(node, tree, tips, assigning, overlap_criterion = 10, coverage_criterion = 0.15):
    leaves1, leaves2 = [child.get_leaves() for child in node.get_children()]
    leaves3 = tips - (set(leaves1) | set(leaves2))
    if not feca2leca(leaves1) or not feca2leca(leaves2) or not feca2leca(leaves3):
        return False
    if assigning:
        species = []
        for leaves in [leaves1, leaves2, leaves3]:
            coverage, copies, repr_species = infer_coverage_redundancy(leaves, representing)
            if coverage < coverage_criterion:
                return False
            species.append(set(repr_species))
        if len(species[0] & species[1]) >= overlap_criterion and len(species[0] & species[2]) >= overlap_criterion and len(species[1] & species[2]) >= overlap_criterion:
            return True
        else:
            return False
    else:
        return True

def duplication_check_rooted(node, assigning, overlap_criterion = 10, coverage_criterion = 0.15):
    daughter1, daughter2 = node.get_children()
    if not duplication_check(daughter1.get_leaves(), daughter2.get_leaves()):
        return False
    if assigning:
        repr_daughter = {}
        repr_daughter[1] = infer_coverage_redundancy(daughter1, representing)
        if repr_daughter[1][0] < coverage_criterion:
            return False
        repr_daughter[2] = infer_coverage_redundancy(daughter2, representing)
        if repr_daughter[2][0] < coverage_criterion:
            return False
        overlap = len(set(repr_daughter[1][2]) & set(repr_daughter[2][2]))
        if overlap < overlap_criterion:
            return False
        else:
            return overlap
    else:
        return True
    
def annotate_and_reroot_euk_only(tree, assigning, overlap_criterion = 10, coverage_criterion = 0.15):
    # Root tree on mid of longest possible distance between the LECAs
    for node in tree.traverse("preorder"):
        node.add_features(identity = "?")
    tips = set(tree.get_leaves())
    dup_counter=0
    for node in tree.iter_descendants("preorder"): ##don't visit root: only two directions (only descendants)
        if not node.is_leaf():
            if duplication_check_unrooted(node, tree, tips, assigning):
                # Duplications only called if there are at least 2 duplications, otherwise no internal node fulfilling this check
                dup_counter += 1
                node.add_features(identity = "duplication", name = "D"+str(dup_counter))
    if dup_counter > 0: # So, at least 2 duplications in the tree
        duplications = tree.search_nodes(identity = 'duplication')
        duplication_paths = {}
        for dup in duplications:
            if dup.get_children()[0] in duplications or dup.get_children()[1] in duplications:
                continue
            name = dup.name
            duplication_paths[name] = [dup.name]
            tree_copy = tree.copy(method = "cpickle")
            search_node = tree_copy&name
            while search_node:
                parent = search_node.up
                if parent.identity == "duplication":
                    duplication_paths[name].append(parent.name)
                    search_node = parent
                elif parent.is_root():
                    duplication_paths[name].append("root")
                    del tree_copy
                    break
                else:
                    del tree_copy
                    break
        convergence_point = None
        # Check if there is a continuous path
        for i, dupl in enumerate(duplication_paths):
            if i == 0:
                convergence_point = set(duplication_paths[dupl])
            else:
                convergence_point = convergence_point & set(duplication_paths[dupl])
        if len(convergence_point) == 0:
            print("Tree contains more than one duplication path with %d duplications in total" % dup_counter, file = sys.stderr)
            sys.exit()
        duplications = tree.search_nodes(identity = "duplication")
        lecas = []
        leca_counter = 0
        for d in duplications:
            daughters = d.get_children()
            if d.up.is_root():
                daughters.append(d.get_sisters()[0])
            else:
                daughters.append(d.up)
            for daughter in daughters:
                if daughter not in duplications:
                    leca_counter += 1
                    name = "OG" + str(leca_counter)
                    daughter.add_features(identity = "LECA", name = name)
                    lecas.append(name)
        ##calculate all distances between LECAs, find the longest distance, root in the middle of this distance
        ##root = duplication
        tree.set_outgroup(duplications[0]) ##tmp rooting on duplication node
        longest_distance = 0
        longest_distance_node1 = ""
        longest_distance_node2 = ""
        for x in range(0, len(lecas)-1):
            for y in range(1, len(lecas)):
                if x < y:
                    distance = tree.get_distance(lecas[x], lecas[y])
                    if distance > longest_distance:
                        longest_distance = distance
                        longest_distance_node1 = lecas[x]
                        longest_distance_node2 = lecas[y]
        mid_longest_distance = 0.5 * longest_distance
        ##find path between node1 and node2 via their last common ancestor
        node1 = tree&longest_distance_node1
        node2 = tree&longest_distance_node2
        lca = tree.get_common_ancestor(node1, node2)
        path = collections.OrderedDict() # Ordered dictionary in which the distance between nodes are stored
        while node1:
            if node1 != lca:
                path[node1] = tree.get_distance(node1, node1.up)
                node1 = node1.up
            else:
                break
        path2 = collections.OrderedDict()
        while node2:
            if node2 != lca:
                path2[node2.up] = tree.get_distance(node2, node2.up)
                node2 = node2.up
            else:
                break
        path2_reverse = collections.OrderedDict(reversed(path2.items()))
        path.update(path2_reverse)
        ##make a list of the path
        node_length_path = []
        for key, value in path.items():
            node_length_path.extend([key, value])
        node_length_path.append(tree&longest_distance_node2) ##complete the path list
        ##find root position
        nodes_on_path = node_length_path[0::2] ##take nodes only
        parent_node_root = ""
        child_node_root = ""
        parent_node_distance = 0
        child_node_distance = 0
        dist_sum = 0
        for x in range(len(nodes_on_path)):
            n = nodes_on_path[x]
            dist = node_length_path[x * 2 + 1]
            dist_sum += dist
            if dist_sum > mid_longest_distance:
                n1 = nodes_on_path[x + 1]
                if n1 == n.up:
                    child_node_root = n
                    parent_node_root = n1
                    parent_node_distance = dist_sum - mid_longest_distance
                else:
                    child_node_root = n1
                    parent_node_root = n
                    parent_node_distance = mid_longest_distance - (dist_sum - dist)
                child_node_distance = dist - parent_node_distance
                break
        ##set root position
        child_node_root_detached = child_node_root.detach()
        child_node_root_detached.add_features(dist = child_node_distance) 
        parent_node_root.add_child(name = "R", dist = parent_node_distance)
        R = tree&"R"
        R.add_child(name = "O", dist = 0) ##outgroup
        R.add_child(child_node_root_detached)
        tree.set_outgroup("O")
        tree = R.detach()
        R.add_features(identity = "duplication", name = "D" + str(dup_counter + 1))
    else:
        ## Either no duplications or 1 duplication
        ## Try rooting on each internal node
        ## Would the root be a duplication node? 
        ## What is the number of species overlap? 
        ## Define LECAs for maximal species overlap
        tmp_outgroup = ""
        tmp_overlap = 0
        tree_copy = tree.copy(method = 'cpickle') ##keeps original rooting position and lengths
        default_outgroup = tree.get_children()[0]
        tree.set_outgroup(default_outgroup)
        for node in tree.iter_descendants("preorder"):
            if node.is_leaf():
                continue
            tree.set_outgroup(node)
            check = duplication_check_rooted(tree, assigning) # If assigning and fulfilling criteria length overlap, else True/False
            if check:
                if assigning:
                    if check > tmp_overlap:
                        tmp_outgroup = node
                        tmp_overlap = check
                else:
                    if tmp_outgroup == '':
                        tmp_outgroup = node
                    else:
                        print('Multiple rootings possible...', file = sys.stderr)
                        sys.exit()
            ##reset the root
            ##branch lengths to root likely will be altered, but that doesn't matter for now. 
            tree.set_outgroup(default_outgroup)
        if tmp_outgroup != "":
            outgroup = tmp_outgroup
            tree.set_outgroup(outgroup)
            sister_outgroup = outgroup.get_sisters()[0]
            outgroup.add_features(identity = "LECA", name = "OG1")
            sister_outgroup.add_features(identity = "LECA", name = "OG2")
            root = tree.get_tree_root()
            root.add_features(identity = "duplication", name = "D1")
        ##no duplications found: root = LECA
        else:
            tree = tree_copy ##take back the original tree
            tree.add_features(identity = "LECA", name = 'OG1') # For ebl would be nice to have a good root here as well
    return tree

#-----------------------------------------------------------------------------------------------------------------------------------------------

# Parse arguments
optlist, args = getopt.getopt(sys.argv[1:], 't:p:o:b:efirh')
opts = {}
euk_only = False; farthest = False; filtering = False; bbh_annotation = False; assigning = False # Default values
if len(args) != 0:
    print('Error: not all arguments recognised', file = sys.stderr); usage()
for k, v in optlist:
    if k == '-h':
        usage()
    else:
        opts[k] = v
if '-t' not in opts:
    print('Error: specify a tree file (-t)', file = sys.stderr); usage()
if '-p' in opts:
    prefix = opts['-p']
else:
    prefix = opts['-t'][:opts['-t'].find('.')]
if '-o' in opts:
    outdir = opts['-o']
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    outdir_prefix = outdir + '/' + prefix
else:
    outdir_prefix = prefix
if '-b' in opts:
    bbh_annotation = True
    try:
        bbh_pairs_file = open(opts['-b'])
    except IOError:
        print('Error:', opts['-b'], 'cannot be opened', file = sys.stderr); sys.exit()
    print('Warning: OG annotation based on BBHs is experimental; BBH triangles not included yet...', file = sys.stderr)
    og_annotation = open(outdir_prefix + '_og_annotation.log', 'w')
if '-e' in opts:
    euk_only = True
if '-f' in opts:
    farthest = True
    print('Warning: only farthest not fully implemented yet', file = sys.stderr)
if '-i' in opts:
    filtering=True
    if euk_only:
        print('Error: filtering of prokaryotes not possible for eukaryote-only trees', file = sys.stderr); usage()
if '-r' in opts:
    assigning = True

# Open tree, get supergroups and annotate leaves
tree = PhyloTree(opts['-t'], sp_naming_function = None)
supergroups2, supergroups5 = get_supergroups()
annotate_prokaryotic_eukaryotic_leaves(tree, euk_only)

# Assign all original sequences to representing sequence
if assigning:
    if euk_only:
        all_seqs_file = open('/hosts/linuxhome/memory/julian2/pfam_hmm/domains_extracted/combined/t1/euk_only/' + prefix + '_seqids.list')
    else:
        all_seqs_file = open('/hosts/linuxhome/memory/julian2/pfam_hmm/domains_extracted/combined/t1/' + prefix + '_seqids.list')
    all_seqs = [line.rstrip() for line in all_seqs_file]
    all_seqs_file.close()
    representing, human_represent = assign_all_seqs(tree, all_seqs, euk_only)
    human_seq_info = seq_info(['HSAP'])['HSAP']

# Filter interspersing prokaryotes
if filtering:
    if len(tree.search_nodes(prok_euk = 'Prokaryote')) > 1: # At least 2 prokaryotic sequences in the tree
        to_prune = True
        while to_prune:
            to_prune = prok_filter(tree)
        if len(tree.search_nodes(prok_euk = 'Prokaryote')) == 0: # All prokaryotic sequences filtered
            print('All prokaryotic sequences removed!', file = sys.stderr)
            euk_only = True

# Parse BBH pairs
if bbh_annotation:
    bbh_pairs = []
    for line in bbh_pairs_file:
        line = line.rstrip()
        seqs = line.split(' ')
        if len(seqs) != 2:
            bbh_pairs_file.close()
            print('Error: BBH pairs file does not consist of a pair per line separated by a space', file = sys.stderr); usage()
        bbh_pairs.append((seqs[0], seqs[1]))
    bbh_pairs_file.close()

# Open output files and print headers
duplication_lengths_out = open(outdir_prefix + '_duplication_lengths.tsv', 'w')
lecas_out = open(outdir_prefix + '_lecas.tsv', 'w')
if not euk_only:
    ancestry_out = open(outdir_prefix + '_ancestry.tsv', 'w')
    print('Pfam\tFECA\tFECA support\tLECAs\tUnknowns\tSister1 name\tSister1 support\tSister2 name\tSister2 support\tAncestry\tSupport', file = ancestry_out)
    branch_lengths_out = open(outdir_prefix + '_branch_lengths.tsv', 'w')
    print('Pfam\tFECA\tAncestry\tLECAs\tProkaryotic sister branch length (med)\tRaw stem lengths (med)\tStem lengths (med)\tEukaryotic branch lengths (med(med))', file = branch_lengths_out)
    non_feca_out = open(outdir_prefix + '_non_feca.tsv', 'w')
if assigning:
    human_seqs = [seq for seq in all_seqs if 'HSAP' in seq]
    if euk_only:
        print('Pfam\tDuplication\tSpecies overlap\tOGs\tRaw duplication lengths (med)\tDuplication lengths (med)', file = duplication_lengths_out)
        print('Pfam\tLECA\tCoverage\tCopy number\tHuman seqs\tHuman name\tSeqs', file = lecas_out)
    else:
        print('Pfam\tFECA\tAncestry\tDuplication\tSpecies overlap\tOGs\tRaw duplication lengths (med)\tDuplication lengths (med)', file = duplication_lengths_out)
        print('Pfam\tFECA\tAncestry\tLECA\tCoverage\tCopy number\tHuman seqs\tHuman name\tSeqs', file = lecas_out)
        print('Pfam\tDonor taxon\tAncestry\tSupport\tCoverage\tSpecies\tSequence IDs', file = non_feca_out)
else:
    if euk_only:
        print('Pfam\tDuplication\tRaw duplication lengths (med)\tDuplication lengths (med)', file = duplication_lengths_out)
        print('Pfam\tLECA\tSeqs', file = lecas_out)
    else:
        print('Pfam\tFECA\tAncestry\tDuplication\tRaw duplication lengths (med)\tDuplication lengths (med)', file = duplication_lengths_out)
        print('Pfam\tFECA\tAncestry\tLECA\tSeqs', file = lecas_out)
        print('Pfam\tDonor taxon\tAncestry\tSupport\tSpecies\tSequence IDs', file = non_feca_out)

# Perform analysis
if euk_only:
    # First check if the tree fulfills the FECA-2-LECA criterion (always the case for 2 supergroups trees)
    if not feca2leca(tree.get_leaves()):
        print('No LECAs in this tree', file = sys.stderr)
        sys.exit()
    # Midpoint rooting to have a root to work with
    root = tree.get_midpoint_outgroup()
    tree.set_outgroup(root)
    # Call duplications in unrooted mode and reroot tree
    tree = annotate_and_reroot_euk_only(tree, assigning)
    if assigning:
        if not annotate_overlap_all_assigned(tree): # Annotate nodes and get relevant coverage and redundancy values
            print('No confident LECA in this tree: coverage =', tree.coverage, file = sys.stderr)
            sys.exit()
    for leca in tree.iter_search_nodes(identity = 'LECA'):
        seqs = [leaf.name for leaf in leca.iter_leaves()]
        if assigning:
            leca_human = {}
            found_human = []
            for human in human_seqs:
                if human in human_represent:
                    if human_represent[human] in seqs:
                        if 'domain' in human:
                            human_seq = human[:human.find('_')]
                        else:
                            human_seq = human
                        leca_human[human] = human_seq_info[human_seq][0]
                        found_human.append(human)
            for human in found_human:
                human_seqs.remove(human)
            print(prefix, leca.name, str(leca.coverage), str(leca.redundancy), ','.join(list(leca_human.keys())), ','.join(list(leca_human.values())), ','.join(seqs), sep = '\t', file = lecas_out)
        else:
            print(prefix, leca.name, ','.join(seqs), sep = '\t', file = lecas_out)
    # Calculate duplication lengths
    for i, duplication in enumerate(tree.iter_search_nodes(identity = 'duplication')): # What to do with unknowns?!
        dupl_lecas = duplication.search_nodes(identity = 'LECA')
        rdl, dl = calculate_median_duplication_length(tree, duplication, dupl_lecas) # Tree --> euk_clade
        if assigning:
            children = duplication.get_children()
            ogs1 = [leca.name for leca in children[0].iter_search_nodes(identity = 'LECA')]
            ogs2 = [leca.name for leca in children[1].iter_search_nodes(identity = 'LECA')]
            print(prefix, duplication.name, duplication.overlap, ','.join(ogs1) + ' - ' + ','.join(ogs2), rdl, dl, sep = '\t', file = duplication_lengths_out)
        else:
            print(prefix, duplication.name, rdl, dl, sep = '\t', file = duplication_lengths_out)
    no_lecas = len(tree.search_nodes(identity = 'LECA'))
    #no_unknowns = len(tree.search_nodes(identity = 'unknown'))
    #print('Pfam\tLECAs\tUnknowns')
    print('Pfam\tLECAs')
    #print(prefix, no_lecas, no_unknowns, sep = '\t')
    print(prefix, no_lecas, sep = '\t')
else:
    euk_clades = get_euk_clades(tree)
    feca_clades = []
    non_feca_clades = []
    feca_count = 0
    for euk_clade in euk_clades:
        if feca2leca(euk_clade.get_leaves()):
            feca_clades.append(euk_clade) # To prevent reuse of same euk clades in the loop
        else:
            non_feca_clades.append(euk_clade)
    if assigning:
        failed_fecas = []
        for euk_clade in feca_clades:
            feca_confidence = annotate_overlap_all_assigned(euk_clade, overlap_criterion = 10, coverage_criterion = 0.15)
            if not feca_confidence:
                failed_fecas.append(euk_clade)
                non_feca_clades.append(euk_clade)
        for failed_feca in failed_fecas:
            feca_clades.remove(failed_feca)
    else:
        annotate_eukaryotic_nodes(feca_clades)
    for euk_clade in feca_clades:
        feca_count += 1
        feca_no = 'FECA' + str(feca_count)
        if bbh_annotation:
            print(feca_no, file = og_annotation)
            reannotate_bbhs(euk_clade, bbh_pairs)
        euk_clade.add_features(feca = True, feca_no = feca_no)
        feca_support = euk_clade.support
        tree_copy = tree.copy(method = 'cpickle') # Copy to prevent errors due to rerooting on sequence of other FECA
        feca_clade = tree_copy.search_nodes(feca_no = feca_no)[0]
        lecas = feca_clade.search_nodes(identity = 'LECA')
        unknowns = feca_clade.search_nodes(identity = 'unknown')
        support = None
        ancestry = ''
        if len(tree_copy.search_nodes(prok_euk = 'Prokaryote')) == 1:
            prok = tree_copy.children[0] # As the tree is rooted on this single prokaryote
            ancestry = classify_sister(prok.taxid)
            lcas = [(prok.taxid, prok.sp, 'NA'),('NA', 'NA', 'NA')]
            support = 'NA'
            feca_support = 'NA'
            print('Warning: one prokaryotic sequence for %s, so branch length analysis may not be accurate' % prefix, file = sys.stderr)
        else:
            lcas = get_prokaryotic_sister(feca_clade, tree_copy, farthest)
            ancestries = []
            for lca in lcas:
                pot_ancestry = classify_sister(lca[0])
                ancestries.append(pot_ancestry)
            if ancestries[0] == ancestries[1]:
                ancestry = ancestries[0]
                if 'NA' in (lcas[0][2], lcas[1][2]):
                    support = 'NA'
                else:
                    support = max(lcas[0][2], lcas[1][2])
                farthest_leaf = reroot(feca_clade, tree_copy) # Root on the farthest leaf
                print('Warning: for %s in %s not enough information to choose a good outgroup and therefore rooted on farthest leaf %s' % (feca_no, prefix, farthest_leaf), file = sys.stderr)
            else:
                order = ['alphaproteobacterial', 'Asgardian', 'aerobic proteobacterial', 'Asgardian+TACK', 'betaproteobacterial', 'gammaproteobacterial', 'TACK']
                for ancestry in order:
                    if ancestry in ancestries:
                        index = ancestries.index(ancestry)
                        support = lcas[index][2]
                        if index == 1:
                            root = feca_clade.get_sisters()[0] # Root on the old sister
                            tree_copy.set_outgroup(root)
                        break
                else:
                    old_sister_leaves = [leaf.name for leaf in feca_clade.get_sisters()[0]]
                    farthest_leaf = reroot(feca_clade, tree_copy) # Root on the farthest leaf
                    print('Warning: for %s in %s not enough information to choose a good outgroup and therefore rooted on farthest leaf %s' % (feca_no, prefix, farthest_leaf), file = sys.stderr)
                    if farthest_leaf in old_sister_leaves: # To check in which of the two possible sisters the farthest leaf is --> other is the sister group
                        index = 1
                    else:
                        index = 0
                    ancestry = ancestries[index]
                    support = lcas[index][2]
        print(prefix, feca_no, feca_support, len(lecas), len(unknowns), lcas[0][1], lcas[0][2], lcas[1][1], lcas[1][2], ancestry, support, sep = '\t', file = ancestry_out)

        # Calculate branch lengths
        sister = feca_clade.get_sisters()[0]
        if len(sister.search_nodes(identity = 'LECA')) > 0:
               print('Warning: for %s in %s there is another FECA in the sister group and therefor the branch lengths may not be accurate' % (feca_no, prefix), file = sys.stderr)
               # To do --> merge these FECAs
        lepca = feca_clade.up
        psbl, rsl, sl, ebl = calculate_branch_lengths(tree_copy, lecas, lepca, sister)
        print(prefix, feca_no, ancestry, len(lecas), psbl, rsl, sl, ebl, sep = '\t', file = branch_lengths_out)

        for i, leca in enumerate(feca_clade.iter_search_nodes(identity = 'LECA')):
            leca_id = 'OG' + str(feca_count) + '.' + str(i+1)
            leca.add_features(og = leca_id)
            seqs = [leaf.name for leaf in leca.iter_leaves()]
            if assigning:
                leca_human = {}
                found_human = []
                for human in human_seqs:
                    if human in human_represent:
                        if human_represent[human] in seqs:
                            if 'domain' in human:
                                human_seq = human[:human.find('_')]
                            else:
                                human_seq = human
                            leca_human[human] = human_seq_info[human_seq][0]
                            found_human.append(human)
                for human in found_human:
                    human_seqs.remove(human)
                print(prefix, feca_no, ancestry, leca_id, str(leca.coverage), str(leca.redundancy), ','.join(list(leca_human.keys())), ','.join(list(leca_human.values())), ','.join(seqs), sep = '\t', file = lecas_out)
            else:
                print(prefix, feca_no, ancestry, leca_id, ','.join(seqs), sep = '\t', file = lecas_out)
        # Calculate duplication lengths
        for i, duplication in enumerate(feca_clade.iter_search_nodes(identity = 'duplication')): # What to do with unknowns?!
            dupl_lecas = duplication.search_nodes(identity = 'LECA')
            rdl, dl = calculate_median_duplication_length(tree_copy, duplication, dupl_lecas) # Tree --> euk_clade
            dupl_id = 'D' + str(feca_count) + '.' + str(i+1)
            if assigning:
                children = duplication.get_children()
                ogs1 = [leca.og for leca in children[0].iter_search_nodes(identity = 'LECA')]
                ogs2 = [leca.og for leca in children[1].iter_search_nodes(identity = 'LECA')]
                print(prefix, feca_no, ancestry, dupl_id, duplication.overlap, ','.join(ogs1) + ' - ' + ','.join(ogs2), rdl, dl, sep = '\t', file = duplication_lengths_out)
            else:
                print(prefix, feca_no, ancestry, dupl_id, rdl, dl, sep = '\t', file = duplication_lengths_out)
    for count,euk_clade in enumerate(non_feca_clades):
        euk_clade.add_features(non_feca_no = count)
        tree_copy = tree.copy(method = 'cpickle')
        non_feca_clade = tree_copy.search_nodes(non_feca_no = count)[0]
        lca, lca_name, support = get_non_feca_sister(non_feca_clade, tree_copy)
        ancestry = classify_sister(lca)
        species = []
        sequences = []
        if assigning:
            supergroup_counts = {'Obazoa':0, 'Amoebozoa':0, 'RASH':0, 'Archaeplastida':0,'Excavata':0}
            repr_species = []
        for leaf in non_feca_clade:
            sp = leaf.taxid
            seqid = leaf.name
            if assigning:
                if seqid in representing:
                    repr_seqs = representing[seqid] + [seqid]
                else:
                    repr_seqs = [seqid]
                for seq in repr_seqs:
                    spec = seq[0:4]
                    if spec not in repr_species:
                        supergroup_counts[supergroups5[spec]] += 1
                        repr_species.append(spec)
            species.append(sp)
            sequences.append(seqid)
        if assigning:
            coverages = []
            number_species = {'Obazoa':123, 'Amoebozoa':6, 'RASH':43, 'Archaeplastida':24,'Excavata':13}
            for group in supergroup_counts:
                coverage = supergroup_counts[group] / number_species[group]
                coverages.append(coverage)
            coverage = sum(coverages) / 5
            print(prefix, lca_name, ancestry, support, coverage, ','.join(species), ','.join(sequences), sep = '\t', file = non_feca_out)
        else:
            print(prefix, lca_name, ancestry, support, ','.join(species), ','.join(sequences), sep = '\t', file = non_feca_out)
    no_lecas = len(tree.search_nodes(identity = 'LECA'))
    no_unknowns = len(tree.search_nodes(identity = 'unknown'))
    print('Pfam\tFECAs\tLECAs\tUnknowns\tNon-FECAs')
    print(prefix, feca_count, no_lecas, no_unknowns, len(non_feca_clades), sep = '\t')
    if assigning:
        for human_seq in human_seqs: # Ones not assigned
            if 'domain' in human_seq:
                human_seqid = human_seq[:human_seq.find('_')]
            else:
                human_seqid = human_seq
            info = human_seq_info[human_seqid][0]
            if human_seq in human_represent:
                print(human_seq, info, human_represent[human_seq], sep = '\t')
            else:
                print(human_seq, info, 'unassigned', sep = '\t')

# Close output files
if bbh_annotation:
    og_annotation.close()
if not euk_only:
    ancestry_out.close()
    branch_lengths_out.close()
    non_feca_out.close()
duplication_lengths_out.close()
lecas_out.close()

# To do:
# If LCA of sister group is a higher level than phylum --> pick majority phylum / proteobacterial class or defined superphylum / superclass
