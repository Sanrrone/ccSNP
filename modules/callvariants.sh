set -e

function callvariants {
	
	echo "$(tput setaf 2)ccSNP: call variants step start$(tput sgr0)"

	SCALLER=$1
	BAMFILES=$2
	REFERENCE=$3
	OUTPUT=$4
	PLOIDY=$5
	QUALITY=$6



	cd $(dirname $REFERENCE)
		refname=$(basename $REFERENCE | awk -F'.' '{for(i=1;i<NF-1;i++){printf("%s.",$i)}printf("%s\n",$i)}')
		samtools dict $REFERENCE -o ${refname}.dict
	cd ..

	mkdir -p ${OUTPUT}_3-callvariant
	cd ${OUTPUT}_3-callvariant



	for scaller in $(echo $SCALLER | tr "," "\n")
	do
		for bam in $(echo "$BAMFILES" | tr ',' ' ')
		do

			case $scaller in
			    gatk|GATK)
					echo "$(tput setaf 2)GATK: running gatk $bam$(tput sgr0)"
					#GCYCLES is defined by ccSNP parameter
					gatkPipeline $REFERENCE $bam $PLOIDY $QUALITY $GCYCLES
			    ;;
			    samtools|SAMTOOLS|SAMtools)
			    ;;
			    freebayes|FREEBAYES|fb|FB)
			    ;;
			esac
		done
	done

	cd ..
	export VCFFILES=$(ls $(pwd)/${OUTPUT}_3-callvariant/*.vcf)
	echo "$(tput setaf 2)ccSNP: call variants step done$(tput sgr0)"
}

#create trusted vcf




#make consensus sequence

#for vcf in $(ls *.vcf)
#do
# sampleName=$(basename $vcf | sed 's/.vcf//g')
# bcftools view $vcf -Oz -o ${sampleName}.bcf
# bcftools index --threads "$(nproc)" ${sampleName}.bcf
# REF=$(echo $sampleName | awk -F"_" '{gsub("ref","",$1);print $1".fasta"}')
# echo "consensus on $REF - $sampleName"
# bcftools consensus -f ../../2-ref/$REF ${sampleName}.bcf > consensus_${sampleName}.fna
# sed -i "s/>.*/>${sampleName}/" consensus_${sampleName}.fna
#done

#statistics
#for vcf in $(ls *realigned.vcf); do grep -v "#" $vcf |grep INDEL |wc -l; done #for indels
#for vcf in $(ls *realigned.vcf); do grep -v "#" $vcf |grep -v INDEL |wc -l; done #for snp only
