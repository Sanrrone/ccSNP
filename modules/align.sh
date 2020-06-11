function align {
set -e

ALIGNER=$1
READS=$2
REFERENCE=$3
OUTPUT=$4

REFIDX=$(basename $REFERENCE | awk '{n=split($1,a,".");for(i=1;i<n;i++){if(i<n-1){printf("%s.",a[i])}else{printf("%s",a[i])}};printf "\n"}')

	echo "$(tput setaf 2)ccSNP: align step start$(tput sgr0)"

	mkdir -p reference
	cd reference
		cp $REFERENCE .
		REFERENCE=$(echo $(pwd)/$(basename $REFERENCE))
		samtools faidx $REFERENCE
	cd ..
	mkdir -p ${OUTPUT}_1-map

	for aligner in $(echo $ALIGNER | tr "," "\n")
	do
			###################################	Indexing zone
			cd reference
				case $aligner in
					bwa|BWA)
						echo "$(tput setaf 2)ccSNP: indexing $REFERENCE using BWA$(tput sgr0)"
						bwa index $REFERENCE
						#REFIDX=$(echo "$(pwd)/$REFIDX")
				    ;;
				    smalt|SMALT)
						echo "$(tput setaf 2)ccSNP: indexing $REFERENCE using SMALT$(tput sgr0)"
						refname=$(basename $REFERENCE | awk -F'.' '{for(i=1;i<NF-1;i++){printf("%s.",$i)}printf("%s\n",$i)}')
						smalt index ${refname} $REFERENCE
				    ;;
				    novoalign|NOVOALIGN)
						novoindex $REFERENCE.idx $REFERENCE
				    ;;
				esac
			cd ..

			################################### Mapping zone
			echo $READS | awk -F';' '{n=split($1,r1,",");split($2,r2,",");for(i=1;i<=n;i++){print r1[i]"\t"r2[i]}}' | while read R1 R2
			do
				sampleName=$(recognizeSampleName $(basename $R1) $(basename $R2))
				reads="$(readlink -f ../$R1) $(readlink -f ../$R2)"
					cd ${OUTPUT}_1-map
					    case $aligner in
					        bwa|BWA)
								echo "$(tput setaf 2)BWA: running bwa mem on $reads$(tput sgr0)"
								bwa mem  -Y -M -R "@RG\tID:${sampleName}\tSM:${sampleName}\tPL:ILLUMINA\tPG:bwa" -t $(nproc) $REFERENCE $reads | samtools view -b -o ${sampleName}.bam -@ $(nproc)
					        ;;
					        smalt|SMALT)
 								smalt map -f sam -o ${sampleName}.smalt.bam -q 10 -n $(nproc) $REFERENCE $reads
					        ;;
					        novoalign|NOVOALIGN)
								novoalign -f $REFERENCE -d $ref -c $(nproc) -o SAM | samtools view -b -o ${sampleName}.novoalign.bam -@ $(nproc)
					        ;;
					    esac
				    cd ..
			done
	done


export BAMFILES=$(ls -1 $(pwd)/${OUTPUT}_1-map/*.bam | tr '\n' ',')
echo "$(tput setaf 2)ccSNP: align step done$(tput sgr0)"

}