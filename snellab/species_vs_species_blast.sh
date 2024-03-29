#!/bin/bash

# Original script written in Perl by Leny

function usage() {
    echo -e "Usage: \n\tspecies_vs_species_blast.sh <file.fa> <cpus> <outputdir>\n"
    echo "Fasta file should only contain sequences starting with a four-letter abbreviation"
    echo "Dependencies: select_sequences.py, makeblastdb and blastp"
    exit
}

# Check number of arguments
if [ "$#" -ne 3 ]; then
    usage
fi

fasta=$1
prefix=${fasta%.fa}
cpus=$2
outdir=$3

mkdir -p $outdir

# Get species list
grep '>' $fasta | cut -c 2-5 | sort -u > $outdir/${prefix}_ids
list=$(cat $outdir/${prefix}_ids)

# Create file with sequences per species and make a Blast database
for species in $list; do
    speciesfa=$outdir/${prefix}_$species.fa
    select_sequences.py -q $species $fasta > $speciesfa
    makeblastdb -in $speciesfa -dbtype prot
done

# Run all-versus-all-Blastp, excluding your own species. Only one target sequence is reported
for species1 in $list; do
    echo "BLASTing sequences of $species1..."
    for species2 in $list; do
	if [ $species1 != $species2 ]; then
	    blastp -query $outdir/${prefix}_$species1.fa -db $outdir/${prefix}_$species2.fa -max_target_seqs 1 -outfmt 6 -num_threads $cpus >> $outdir/${prefix}_blastp.txt
	fi
    done
done
