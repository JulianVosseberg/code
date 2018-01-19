#!/bin/bash

# Implementation of the strategy advised by Kalyaanamoorthy et al. (2017)

# State usage
function usage() {
    echo -e "Usage: \n\tdoStandardIQTREE.sh -a <alignment_file> [ -p <prefix> ] [ -m <MF+tree|MF|tree(provide model)> ] [ -s <substitution model> ] [ -t <nuclear|mitochondrial|chloroplast|viral> ] [ -c <threads> ] [ -f ] [ -x ]\n"
    echo "Performs ModelFinder and infers phylogeny using the chosen model."
    echo "Default modes are MF+tree, no subtype, auto option for best number of threads and basename alignment file as prefix. Subtype only relevant for ModelFinder."
    echo "If fast option is chosen (-f):"
    echo -e "\t*No advanced search for RHAS model is performed, including check for higher R values if best model is +R10 (both recommended in Kalyaanamorrthy et al. (2017))."
    echo -e "\t*Fast option for tree inference (~FastTree) is performed."
    echo "If a substitution model is provided (-s), only a RHAS model is searched for, either full search (default) or fast search (-f)."
    echo "If extended option is chosen (-x), mixture models (standard only LG: LG+C10, ..., LG+C60; if substitution model is provided, that +C10..60) are tested as well."
    exit
}

#TODO: test for numerical underflow as error (--> safe option)

# If number of arguments is less than 2, invoke usage function
if [ "$#" -lt "2" ]; then
    usage
    exit
fi

# Default mode and subtype
mode="MF+tree"
se_model=""
subtype=""
threads="AUTO"
fast=false
mixture=false

# State options
while getopts ":a:p:m:s:t:c:fx" opt; do
    case $opt in
        a) aln=${OPTARG};;
        p) prefix=${OPTARG};;
        m) mode=${OPTARG};;
	s) se_model=${OPTARG};;
	t) subtype=${OPTARG};;
	c) threads=${OPTARG};;
        f) fast=true;;
	x) mixture=true;;
        *) usage ;;
    esac
done

# Prefix
if [ "$prefix" = "" ]; then
    prefix=$(basename ${aln%.*})
fi

# Interpret mode
if [ "$mode" = "MF" ]; then
    finder=true
    tree=false
elif [ "$mode" = "MF+tree" ]; then
    finder=true
    tree=true
else
    finder=false
    tree=true
    model=$mode
fi

# Get subtype
if [ "$subtype" != "" ]; then
    msub="-msub $subtype"
else
    msub=""
fi

# Get substitution model
if [ "$se_model" != "" ]; then
    mset=true
else
    mset=false
fi

# Get mixture models
if $mixture; then
    if ! $mset; then
	se_model="LG"
    fi
    madd="-madd C10+$se_model,C20+$se_model,C30+$se_model,C40+$se_model,C50+$se_model,C60+$se_model"
else
    madd=""
fi

# Perform ModelFinder
if $finder; then
    if ! $mset; then
	echo "Performing ModelFinder..."
	iqtree -s $aln -m MF -pre ${prefix}_MF -nt $threads $msub $madd -quiet
	model=$(grep '^Best-fit model according to BIC:' ${prefix}_MF.iqtree | sed 's/Best-fit model according to BIC: //')
	if [ "$model" = "" ]; then
	    echo "Error: Best-fit model not detected"
	    exit
	fi
	echo -e "Best-fit model is $model\n"
    fi
    if ! $fast; then
       echo "Performing advanced search for the optimal model of rate heterogeneity across sites..."
       se_model=${model%%+*}
       if [[ "$se_model" = "C"[1-6]"0" ]]; then
	   echo "Skipped: no advanced search for RHAS in case of a mixture model"
       else 
	   iqtree -s $aln -pre ${prefix}_MF_$se_model -mset $se_model -nt $threads -m MF -mtree -quiet
	   model=$(grep '^Best-fit model according to BIC:' ${prefix}_MF_$se_model.iqtree | sed 's/Best-fit model according to BIC: //')
	   if [ "$model" = "" ]; then
	       echo "Error: Best-fit model not detected"; exit
	   fi
       fi
       echo -e "Best-fit model is $model\n"
       if [[ "$model" = *"+R10" ]]; then
	  se_model=${model%%+*}
	  echo "Best-fit model includes the R10 model of RHAS, also checking $se_model models with higher R values..."
	  iqtree -s $aln -pre ${prefix}_MF_${se_model}_R8-20 -mset $se_model -nt $threads -m MF -mtree -cmin 8 -cmax 20 -quiet
	  model=$(grep '^Best-fit model according to BIC:' ${prefix}_MF_${se_model}_R8-20.iqtree | sed 's/Best-fit model according to BIC: //')
	  if [ "$model" = "" ]; then
	    echo "Error: Best-fit model not detected"; exit
	  fi
	  echo -e "Best-fit model is $model\n"
       fi
    fi
    if $fast && $mset; then
	echo "Performing fast search for the optimal $se_model model..."
	iqtree -s $aln -pre ${prefix}_MF_$se_model -mset $se_model -nt $threads -m MF -quiet $madd
	model=$(grep '^Best-fit model according to BIC:' ${prefix}_MF_$se_model.iqtree | sed 's/Best-fit model according to BIC: //')
	if [ "$model" = "" ]; then
	    echo "Error: Best-fit model not detected"; exit
	fi
	echo -e "Best-fit model is $model\n"
    fi
fi

# Infer phylogenetic tree
if $tree; then
    if $fast; then
	echo "Inferring fast tree with $model model..."
	iqtree -s $aln -m $model -pre $prefix -nt $threads -alrt 1000 -quiet -fast
    else
	echo "Inferring tree with $model model..."
	iqtree -s $aln -m $model -pre $prefix -nt $threads -alrt 1000 -bb 1000 -quiet
    fi
    echo "Done!"
fi
