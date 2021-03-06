#!/usr/bin/env python3

# Functions for using FECA-to-LECA ScrollSaw trees

from numpy import median

def assign_all_seqs(tree_seqs, all_seqs, euk_only, prefix, blast_path = '/home/julian/julian2/pfam_hmm/improved_pipeline/euk_fasta/blast_output'):
    """Assigns all original sequences to their best representing tree sequence based on the BLAST results and keeps track of human sequences"""
    other_seqs = set(all_seqs) - set(tree_seqs)
    human_represent = {}
    for seq in tree_seqs:
        if 'HSAP' in seq:
            human_represent[seq] = seq # Human seq represented by itself
    with open(f'{blast_path}/{prefix}/{prefix}_blastp.txt') as blast:
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
    # To correctly assign in-paralogues: also check the own vs. own BLAST
    with open(f'{blast_path}/{prefix}/{prefix}_blastp_own.txt') as blast_own:
        for line in blast_own:
            line = line.rstrip()
            fields = line.split('\t')
            query, hit = fields[0:2]
            if query == hit or query in tree_seqs or hit not in tree_seqs:
                continue
            bit_score = float(fields[11])
            if query not in best_blast_hits:
                best_blast_hits[query] = (hit, bit_score)
            else:
                best_score = best_blast_hits[query][1]
                if bit_score > best_score:
                    best_blast_hits[query] = (hit, bit_score)
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

def infer_coverage_redundancy(node, representing, supergroups, tree = True):
    number_species = {}
    for group in supergroups.values():
        number_species[group] = number_species.get(group, 0) + 1
    all_species = []
    species_counts = {}
    supergroup_counts = {group: 0 for group in number_species}
    group_species_counts = {group: [] for group in number_species}
    for leaf in node: # node can also be a list of leaves instead
        if tree:
            name = leaf.name
            sp = leaf.taxid
        else:
            name = leaf
            sp = name[0:4]
        species = [seq[0:4] for seq in representing.get(name,'')]
        species.append(sp)
        all_species.extend(species)
    for spec in all_species:
        if spec in species_counts:
            species_counts[spec] += 1
        else:
            supergroup = supergroups[spec]
            supergroup_counts[supergroup] += 1
            species_counts[spec] = 1
    coverages = []
    for group, counts in supergroup_counts.items():
        coverages.append(counts / number_species[group])
    coverage_av = sum(coverages) / len(number_species)
    for spec, counts in species_counts.items():
        group = supergroups[spec]
        group_species_counts[group].append(counts)
    redundancy_med = median([median(counts) for counts in group_species_counts.values() if len(counts) > 0])
    return coverage_av, redundancy_med, list(species_counts.keys())

def get_human_representing(seqs, human_seqs, human_represent, human_seq_info):
    remaining_human_seqs = human_seqs[:]
    repr_human = {}
    for human in human_seqs:
        if human in human_represent:
            if human_represent[human] in seqs:
                if '_' in human:
                    human_seq = human[:human.find('_')]
                else:
                    human_seq = human
                repr_human[human] = human_seq_info[human_seq][0]
                remaining_human_seqs.remove(human)
    return remaining_human_seqs, repr_human
