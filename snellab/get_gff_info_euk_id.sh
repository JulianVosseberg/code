fasta=$1

if [[ -f gff_info.txt ]]; then
    echo "Error: gff_info.txt already exists" >&2
    exit
fi

ids=$(grep '>' $fasta | sed -r 's/>.+_//')
for id in $ids; do
    abbr=$(echo $id | sed -r 's/[0-9]+//')
    origid=$(grep -m 1 $id ~/julian2/snel-clan-genomes/eukarya/proteomes_metadata/$abbr.metadata.txt | cut -f 11)
    if [[ $origid = "ENS"* ]]; then
	origid=${origid%.*}
    fi
    if ! [[ -f ~/julian2/snel-clan-genomes/eukarya/gff/$abbr.gff3.gz || -f ~/julian2/snel-clan-genomes/eukarya/gff/$abbr.gff.gz ]]; then
	echo "No gff file for $abbr" >&2
    else
	echo -e "$abbr\t$id\t$origid" >> gff_info.txt
	pattern="$origid;"
	if [[ $abbr = "ODIO" ]]; then
	    origid=GSOID_T${origid#GSOIDP}
	    origid=rna$(zcat ~/julian2/snel-clan-genomes/eukarya/gff/ODIO.gff.gz | grep $origid | sed -r 's/.+ID=gene//; s/;Name.+//')
	    pattern="$origid;"
	elif [[ $abbr = "BSCH" ]]; then
            origid=$(echo $origid | sed 's/_/.t/')
            pattern="$origid\";"
        elif [[ $abbr = "GINT" || $abbr = "MONO" || $abbr = "PPAC" || $abbr = "SMIN" ]]; then
            pattern="$origid$"
        elif [[ $abbr = "NGAD" ]]; then
            pattern="$origid "
        elif [[ $abbr = "PMIN" ]]; then
            pattern=${origid%-RA}-tr
        fi
	zcat ~/julian2/snel-clan-genomes/eukarya/gff/$abbr.gff*.gz | grep $pattern >> gff_info.txt
    fi
done

# Check non-retrieved information:
## HSAP, CINT, DRER, XTRO, APLA, PSIN, OANA: Ensembl IDs --> exclude decimal at end
## ODIO: GSOIDP# in metadata file --> GSOID_T# in gff --> ID=gene# in same line --> grep rna# lines
## ESIL: also not found manually in gff file; identifiers not corresponding
