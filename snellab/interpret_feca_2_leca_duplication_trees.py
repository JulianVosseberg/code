#!/home/julian/.local/anaconda3/bin/python3

# Load modules
from ete3 import PhyloTree
from ete3 import NCBITaxa
from numpy import median
import re
import sys
import getopt
import os
ncbi = NCBITaxa()

# Usage
def usage():
    print('\n\tUsage: interpret_feca_2_leca_duplication_trees.py -t <tree.contree> [ -p <prefix> ] [ -o <output_dir> ] [ -e ] [ -f ] [ -h ]', file = sys.stderr)
    print('\nThis script identifies FECA-2-LECA duplications, determines the best prokaryotic outgroup (if any), and performs a branch length analysis', file = sys.stderr)
    print('\nOptions:\n\t-t: tree file\n\t-p: prefix for output files (DEFAULT: basename tree file)\n\t-o: directory for output files (DEFAULT: current)\n\t-e: only eukaryotes (DEFAULT: off)\
    \n\t-f: use only farthest leaf for rooting (DEFAULT: off)\n\t-h: help', file = sys.stderr)
    sys.exit()

# Parse supergroups
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

# Distinguish eukaryotic and prokaryotic leaves
def annotate_prokaryotic_eukaryotic_leaves(tree, euk_only):
    if euk_only:
        for leaf in tree:
            taxid = leaf.name[0:4]
            leaf.add_features(supergroup2 = supergroups2[taxid], supergroup5 = supergroups5[taxid])
    else:
        for leaf in tree:
            if re.match('^\d', leaf.name):
                taxid = int(leaf.name[:leaf.name.find('.')])
                lineage = ncbi.get_lineage(taxid)
                leaf.add_features(taxid = taxid, sp = ncbi.get_taxid_translator([taxid])[taxid], prok_euk = 'Prokaryote')
            else:
                taxid = leaf.name[0:4]
                leaf.add_features(taxid = taxid, supergroup2 = supergroups2[taxid], supergroup5 = supergroups5[taxid], prok_euk = 'Eukaryote')

# Count number of monophyletic eukaryotic clades
def get_number_eukaryotic_clades(tree):
    euk_seqs_type = tree.check_monophyly(values = ['Eukaryote'], target_attr = 'prok_euk')[1]
    if euk_seqs_type != 'monophyletic':
        clades = tree.get_monophyletic(values = ['Eukaryote'], target_attr = 'prok_euk')
        count_clades = 0
        for ca in clades:
            count_clades += 1
    else:
        count_clades = 1
    return euk_seqs_type, count_clades

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

# Count number of FECAs
def get_number_fecas(clades):
    feca_count = 0
    for ca in clades:
        if feca2leca(ca.get_leaves()):
            feca_count += 1
    return feca_count

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

#-----------------------------------------------------------------------------------------------------------------------------------------------

# Parse arguments
optlist, args = getopt.getopt(sys.argv[1:], 't:p:o:efh')
opts = {}
euk_only = False; farthest = False # Default values
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
if '-e' in opts:
    euk_only = True
    print('Warning: euk-only not fully implemented yet', file = sys.stderr)
if '-f' in opts:
    farthest = True
    print('Warning: only farthest not fully implemented yet', file = sys.stderr)

# Open tree, get supergroups and annotate leaves
tree = PhyloTree(opts['-t'], sp_naming_function = None)
supergroups2, supergroups5 = get_supergroups()
annotate_prokaryotic_eukaryotic_leaves(tree, euk_only)

# Open output files and print headers
ancestry_out = open(outdir_prefix + '_ancestry.tsv', 'w')
print('Pfam\tFECA\tFECA support\tLECAs\tSister1 name\tSister1 support\tSister2 name\tSister2 support\tAncestry\tSupport', file = ancestry_out)
branch_lengths_out = open(outdir_prefix + '_branch_lengths.tsv', 'w')
print('Pfam\tFECA\tAncestry\tLECAs\tProkaryotic sister branch length (med)\tRaw stem lengths (med)\tStem lengths (med)\tEukaryotic branch lengths (med(med))', file = branch_lengths_out)
duplication_lengths_out = open(outdir_prefix + '_duplication_lengths.tsv', 'w')
print('Pfam\tFECA\tAncestry\tDuplication\tRaw duplication lengths (med)\tDuplication lengths (med)', file = duplication_lengths_out)
non_feca_out = open(outdir_prefix + '_non_feca.tsv', 'w')
print('Pfam\tDonor taxon\tAncestry\tSupport\tSpecies\tSequence IDs', file = non_feca_out)

# Perform analysis
if euk_only:
    pass
else:
    euk_clades = get_euk_clades(tree)
    annotate_eukaryotic_nodes(euk_clades)
    euk_clades = get_euk_clades(tree) # Again? Apparently when iterating they are removed
    feca_clades = []
    non_feca_clades = []
    feca_count = 0
    for euk_clade in euk_clades:
        if feca2leca(euk_clade.get_leaves()):
            feca_clades.append(euk_clade) # To prevent reuse of same euk clades in the loop
        else:
            non_feca_clades.append(euk_clade)
    for euk_clade in feca_clades:
        feca_count += 1
        feca_no = 'FECA' + str(feca_count)
        euk_clade.add_features(feca = True, feca_no = feca_no)
        feca_support = euk_clade.support
        tree_copy = tree.copy(method = 'cpickle') # Copy to prevent errors due to rerooting on sequence of other FECA
        feca_clade = tree_copy.search_nodes(feca_no = feca_no)[0]
        lecas = feca_clade.search_nodes(identity = 'LECA')
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
        print(prefix, feca_no, feca_support, len(lecas), lcas[0][1], lcas[0][2], lcas[1][1], lcas[1][2], ancestry, support, sep = '\t', file = ancestry_out)

        # Calculate branch lengths
        sister = feca_clade.get_sisters()[0]
        if len(sister.search_nodes(identity = 'LECA')) > 0:
               print('Warning: for %s in %s there is another FECA in the sister group and therefor the branch lengths may not be accurate' % (feca_no, prefix), file = sys.stderr)
        lepca = feca_clade.up
        psbl, rsl, sl, ebl = calculate_branch_lengths(tree_copy, lecas, lepca, sister)
        print(prefix, feca_no, ancestry, len(lecas), psbl, rsl, sl, ebl, sep = '\t', file = branch_lengths_out)

        # Calculate duplication lengths
        for i, duplication in enumerate(feca_clade.iter_search_nodes(identity = 'duplication')): # What to do with unknowns?!
            dupl_lecas = duplication.search_nodes(identity = 'LECA')
            rdl, dl = calculate_median_duplication_length(tree_copy, duplication, dupl_lecas) # Tree --> euk_clade
            print(prefix, feca_no, ancestry, 'D'+str(i+1), rdl, dl, sep = '\t', file = duplication_lengths_out)
    for count,euk_clade in enumerate(non_feca_clades):
        euk_clade.add_features(non_feca_no = count)
        tree_copy = tree.copy(method = 'cpickle')
        non_feca_clade = tree_copy.search_nodes(non_feca_no = count)[0]
        lca, lca_name, support = get_non_feca_sister(non_feca_clade, tree_copy)
        ancestry = classify_sister(lca)
        species = []
        sequences = []
        for leaf in non_feca_clade:
            species.append(leaf.taxid)
            sequences.append(leaf.name)
        print(prefix, lca_name, ancestry, support, ','.join(species), ','.join(sequences), sep = '\t', file = non_feca_out)

# Close output files
ancestry_out.close()
branch_lengths_out.close()
duplication_lengths_out.close()
non_feca_out.close()
