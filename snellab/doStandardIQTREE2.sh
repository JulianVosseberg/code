#!/bin/bash

# Implementation of the strategy advised by Kalyaanamoorthy et al. (2017)

# State usage
function usage() {
    echo -e "Usage: \n\tdoStandardIQTREE2.sh -a <alignment_file> [ -p <prefix> ] [ -m <MF+tree|MF|tree(provide model)> ] [ -s <substitution model> ] [ -t <nuclear|mitochondrial|chloroplast|viral> ] [ -c <threads> ] [ -f ] [ -x ] [ -r ]\n"
    echo "Performs ModelFinder and infers phylogeny using the chosen model."
    echo "Default modes are MF+tree, no subtype, auto option for best number of threads and basename alignment file as prefix. Subtype only relevant for ModelFinder."
    echo "If fast option is chosen (-f):"
    echo -e "\t*No advanced search for RHAS model is performed, including check for higher R values if best model is +R10 (both recommended in Kalyaanamorrthy et al. (2017))."
    echo -e "\t*Fast option for tree inference (~FastTree) is performed."
    echo "If a substitution model is provided (-s), only a RHAS model is searched for, either full search (default) or fast search (-[fr])."
    echo "If extended option is chosen (-x), mixture models (standard only LG: LG+C10, ..., LG+C60; if substitution model is provided, that +C10..60; +LG4[MX]) are tested as well."
    echo "If rapid search option is chosen (-r), only standard rapid searches are performed, but including a check for higher R values if best model is R10."
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
advanced=true

# State options
while getopts ":a:p:m:s:t:c:fxr" opt; do
    case $opt in
        a) aln=${OPTARG};;
        p) prefix=${OPTARG};;
        m) mode=${OPTARG};;
	s) se_model=${OPTARG};;
	t) subtype=${OPTARG};;
	c) threads=${OPTARG};;
        f) fast=true;;
	x) mixture=true;;
	r) advanced=false;;
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
    msub="--msub $subtype"
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
    madd="--madd $se_model+C10,$se_model+C20,$se_model+C30,$se_model+C40,$se_model+C50,$se_model+C60,LG4M,LG4X"
else
    madd=""
fi

# Rapid search mode
if $advanced; then
    mtree="--mtree"
else
    mtree=""
fi

# Perform ModelFinder
if $finder; then
    if ! $mset; then
	echo "Performing ModelFinder..."
	iqtree -s $aln -m MF --prefix ${prefix}_MF -T $threads $msub $madd --quiet
	model=$(grep '^Best-fit model according to BIC:' ${prefix}_MF.iqtree | sed 's/Best-fit model according to BIC: //')
	if [ "$model" = "" ]; then
	    echo "Error: Best-fit model not detected"
	    exit
	fi
	echo -e "Best-fit model is $model\n"
    else
	echo "Performing ModelFinder to search for the best-fit $se_model model..."
	iqtree -s $aln -m MF --prefix ${prefix}_MF -T $threads --mset $se_model $madd --quiet
	model=$(grep '^Best-fit model according to BIC:' ${prefix}_MF.iqtree | sed 's/Best-fit model according to BIC: //')
	if [ "$model" = "" ]; then
	    echo "Error: Best-fit model not detected"
	    exit
	fi
	echo -e "Best-fit model is $model\n"
    fi
#    if [[ "$threads" = "AUTO" ]]; then # Not each time auto detect again
#	threads=$(grep -m 1 THREADS ${prefix}_MF.log | sed 's/BEST NUMBER OF THREADS: //')
#	if [ "$threads" = "" ]; then
#	    echo "Error: best number of threads not detected"
#	    exit
#	fi
#    fi
    if ! $fast; then
       if [[ "$model" = *"C"[1-6]"0" ]]; then
	    if $advanced; then
		echo "Performing advanced search for the optimal model of rate heterogeneity across sites..."
	    else
		echo "Performing ModelFinder to search for the best-fit $model model..."
	    fi
	    iqtree -s $aln --prefix ${prefix}_MF_$model --mset $model -T $threads -m MF $mtree --quiet
	    model=$(grep '^Best-fit model according to BIC:' ${prefix}_MF_$model.iqtree | sed 's/Best-fit model according to BIC: //')
       elif $advanced; then
	   echo "Performing advanced search for the optimal model of rate heterogeneity across sites..."
	   if [[ "$model" = "LG4"[MX] ]]; then
	       echo "No RHAS search possible for $model!"
	   else
	       se_model=${model%%+*}
	       iqtree -s $aln --prefix ${prefix}_MF_$se_model --mset $se_model -T $threads -m MF $mtree --quiet
	       model=$(grep '^Best-fit model according to BIC:' ${prefix}_MF_$se_model.iqtree | sed 's/Best-fit model according to BIC: //')
	   fi
       fi
       if [ "$model" = "" ]; then
	   echo "Error: Best-fit model not detected"; exit
       fi
       echo -e "Best-fit model is $model\n"
       if [[ "$model" = *"+R10" ]]; then
	  se_model=${model%%+*}
	  echo "Best-fit model includes the R10 model of RHAS, also checking $se_model models with higher R values..."
	  iqtree -s $aln --prefix ${prefix}_MF_${se_model}_R8-20 --mset $se_model -T $threads -m MF $mtree -cmin 8 -cmax 20 --quiet
	  model=$(grep '^Best-fit model according to BIC:' ${prefix}_MF_${se_model}_R8-20.iqtree | sed 's/Best-fit model according to BIC: //')
	  if [ "$model" = "" ]; then
	    echo "Error: Best-fit model not detected"; exit
	  fi
	  echo -e "Best-fit model is $model\n"
       fi
    fi
fi

# Infer phylogenetic tree
if $tree; then
    if $fast; then
	echo "Inferring fast tree with $model model..."
	iqtree -s $aln -m $model --prefix $prefix -T $threads --alrt 1000 --quiet --fast
    else
	echo "Inferring tree with $model model..."
	iqtree -s $aln -m $model --prefix $prefix -T $threads --alrt 1000 -B 1000 --wbtl --quiet
    fi
    echo "Done!"
fi
