function prepareRef {
set -e

ALIGNER=$1
REFERENCE=$2

	###################################	Indexing zone
	cd reference
		case $ALIGNER in
			bwa|BWA)
				echo "$(tput setaf 2)-- indexing $REFERENCE using BWA$(tput sgr0)"
				bwa index $REFERENCE
				#REFIDX=$(echo "$(pwd)/$REFIDX")
		    ;;
		    smalt|SMALT)
				echo "$(tput setaf 2)-- indexing $REFERENCE using SMALT$(tput sgr0)"
				refname=$(basename $REFERENCE | awk -F'.' '{for(i=1;i<NF-1;i++){printf("%s.",$i)}printf("%s\n",$i)}')
				smalt index ${refname} $REFERENCE
		    ;;
		    novoalign|NOVOALIGN)
				novoindex $REFERENCE.idx $REFERENCE
		    ;;
		esac
	cd ..
}