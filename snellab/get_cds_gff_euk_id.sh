fasta=$1
output=$2
if [[ $output = "" ]]; then
    output=cds_euk_ids.tsv
    if [[ -f $output ]]; then
    echo "Error: $output already exists; please specify another output file name" >&2
    exit
    fi
fi

ids=$(grep '>' $fasta | sed -r 's/>.+_//')
for id in $ids; do
    abbr=$(echo $id | sed -r 's/[0-9]+//')
    origid=$(grep -m 1 $id ~/julian2/snel-clan-genomes/eukarya/proteomes_metadata/$abbr.metadata.txt | cut -f 11)
    if [[ $origid = "ENS"* ]]; then
	origid=${origid%.*} ## HSAP, CINT, DRER, XTRO, APLA, PSIN, OANA: Ensembl IDs --> exclude decimal at end
    fi
    if ! [[ -f ~/julian2/snel-clan-genomes/eukarya/gff/$abbr.gff3.gz || -f ~/julian2/snel-clan-genomes/eukarya/gff/$abbr.gff.gz ]]; then
	echo "No gff file for $abbr" >&2
    else
	pattern="$origid;"
	if [[ $abbr = "ODIO" ]]; then ## ODIO: GSOIDP# in metadata file --> GSOID_T# in gff --> ID=gene# in same line --> grep rna# lines
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
	zcat ~/julian2/snel-clan-genomes/eukarya/gff/$abbr.gff*.gz | grep $pattern | grep -w CDS | cut -f 1,3,4,5,7 | sed -r "s/^.+CDS/$id/" >> $output
    fi
done

# To do:
## ESIL: also not found manually in gff file; identifiers not corresponding
## ADEA/ADAE: implement potential typo
