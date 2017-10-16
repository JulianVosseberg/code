#!/bin/bash

# State usage
function usage() {
    echo -e "Usage: \n\tmake_pthr_alignments.sh -f <file.tsv> [ -c <prefix> ] [ -t <threads> ]\n"
    echo "First step for mapping introns pipeline."
    echo "Based on PTHR hits sequences from Eukarya 4 are selected and aligned."
    echo "A fast tree is inferred to check the assigned sequences."
    echo -e "These trees should be checked by eye, mislabelled and ambiguous sequences should be removed and a new alignment should be made.\n"
    echo "File.tsv should be in the format PTHR#####(.SF#/.orig)\tgene_name(s)."
    echo "Default settings: each PTHR in the file handled separately and not combined with other PTHRs (-c); # threads is 7 (only relevant for MAFFT, for IQTREE auto-detection)"
    exit
}

# If number of arguments is less than 2, invoke usage function
if [ "$#" -lt "2" ]; then
    usage
fi

# Set default mode and values
combined=false # Default: each PTHR in the file separately, not combined
threads="7"

# State options
while getopts ":f:c:t" opt; do
    case $opt in
	f) table=${OPTARG};;
	c) combined=${OPTARG};;
	t) threads=${OPTARG};;
	*) usage ;;
    esac
done

# Directory with for each species a file (ABBR_bestHits) containing the best PTHR hit for each sequence
hits_dir=/hosts/linuxhome/scarab/eva2/pantherOrth/bestHitsPANTHER

# Proteomes file
proteomes=~/julian2/snel-clan-genomes/eukarya/eukarya_proteomes_lt_core.fa

# Make check directory
if [[ ! -d check ]]; then
    mkdir check
else
    echo "Error: Directory 'check' already exists"; exit
fi

# Select sequences
while read line; do
    line=($(echo $line | tr -d '\n' | tr -d '\r'))
    if $combined; then
	selectSequences.pl -i $proteomes -q $(grep ${line[0]} $hits_dir/*_bestHits | sed -r 's/.+://; s/,.+//; s/ADAE/ADEA/' | tr '\n' ',') -t -e | sed "s/>/>${line[1]}_/" >> $combined.fa
    else
	selectSequences.pl -i $proteomes -q $(grep ${line[0]} $hits_dir/*_bestHits | sed -r 's/.+://; s/,.+//; s/ADAE/ADEA/' | tr '\n' ',') -t -e | sed "s/>/>${line[1]}_/" >> ${line[0]%.*}.fa
    fi
done < $table

# Make alignment(s)
if $combined; then
    mafft --localpair --maxiterate 1000 --thread $threads $combined.fa > $combined.aln
else
    for fa in PTHR*.fa; do
	mafft --localpair --maxiterate 1000 --thread $threads $fa > ${fa%fa}aln
    done
fi

# Check for non-orthologues and ambiguous sequences
cd check
if $combined; then
    trimal -in ../$combined.aln -out $combined.gt50.aln -gt 0.5
    doStandardIQTREE.sh -a $combined.gt50.aln -m LG+G4 -f
else
    for aln in ../PTHR*.aln; do
	base_aln=${aln#../}
	trimal -in $aln -out ${base_aln%aln}gt50.aln -gt 0.5
	doStandardIQTREE.sh -a ${base_aln%aln}gt50.aln -m LG+G4 -f
    done
fi
