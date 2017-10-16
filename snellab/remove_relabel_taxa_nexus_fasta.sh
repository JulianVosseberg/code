#!/bin/bash

# State usage
function usage() {
    echo -e "Usage: \n\tremove_relabel_taxa_nexus_fasta.sh <tree.color.nex> <seqs.fa> [ <colours.tsv> ]\n"
    echo "Removes taxa from fasta file based on taxa coloured white in the nexus file and relabels other coloured taxa according to the table provided (colour code (e.g., 00ff00) \t new label)"
    echo "Dependency: selectSequences.pl"
    exit
}

# Check number of arguments
if [ "$#" -lt "2" ]; then
    usage
fi

nexus=$1
fasta=$2
table=$3

# Remove taxa
remov=$(grep '&!color=#ffffff' $nexus | sed 's/\[&!color=#ffffff\]//; s/\t//' | tr '\n', ',')
if [[ $remov = "" ]]; then
    ln -s $fasta ${fasta%.*}_cleaned.fa
else
    selectSequences.pl -i $fasta -q $remov -m remove -e > ${fasta%.*}_cleaned.fa
fi

if [[ $table = "" ]]; then exit; fi

# Relabel taxa
while read line; do
    line=($line)
    colour=${line[0]}
    new_prefix=${line[1]}
    ids=$(grep "&!color=#$colour" $nexus | sed "s/\[&!color=#$colour\]//; s/\t//")
    for id in $ids; do
	new_id=${new_prefix}_${id##*_}
	sed -i "/>/ s/$id/$new_id/" ${fasta%.*}_cleaned.fa
    done
done < $table
	  
