#!/usr/bin/env bash

groups=$1
clan=$2
e_pfams=$3

# Make directory for clan
#mkdir $clan
cd $clan

# Select sequences per FECA per Pfam (in Python script) get_feca_sequence_ids.py
#for pfam in $(grep $clan ../clan_pfams.tsv | cut -f 3 | tr ',' ' '); do
#    ~/proj/timing_dupl/analyses/python_scripts/get_feca_sequence_ids.py $pfam $groups
#done

# Make MSA of FECAs: mafft-linsi
for fa in *fa; do
    ucount=$(grep -v '>' $fa | grep 'U' | head -n 1)
    if [ $ucount ]; then
	mafft-linsi --anysymbol $fa > ${fa%fa}aln
    else
	mafft-linsi $fa > ${fa%fa}aln
    fi
done

# Make HHMs of FECAs: hhmake
for aln in *aln; do
    hhmake -i $aln -o ${aln%aln}hhm -M 50 -name ${aln%.aln}
done

# Make a HHM database: to check
mkdir db; cd db
for hhm in ../*hhm; do
    ln -s $hhm .
done
cd ..
ffindex_build -s ${clan}_FECAs_hhm.ffdata ${clan}_FECAs_hhm.ffindex db/

# Get Pfam HHM queries: get_pfam_hhm.py (stored in pfam_hhm directory)
mkdir euk_only; cd euk_only
get_pfam_hhm.py $(echo $e_pfams | tr ',' ' ')
#$(grep $clan ../../clan_pfams.tsv | cut -f 4 | tr ',' ' ')
cd ..

# Search with Pfam HHMs in FECA HHM database
for hhm in euk_only/*hhm; do
    hhsearch -i $hhm -d ${clan}_FECAs
done
