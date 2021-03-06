set -e

function callvariants {
	
	echo "$(tput setaf 2)ccSNP: call variants step start$(tput sgr0)"

	SCALLER=$1
	BAMFILES=$2
	REFERENCE=$3
	PLOIDY=$4
	QUALITY=$5



	cd $(dirname $REFERENCE)
		refname=$(basename $REFERENCE | awk -F'.' '{for(i=1;i<NF-1;i++){printf("%s.",$i)}printf("%s\n",$i)}')
		samtools dict $REFERENCE -o ${refname}.dict
	cd ..

	mkdir -p 3-callvariant
	cd 3-callvariant



	for bam in $(echo "$BAMFILES" | tr ',' ' ')
	do
		sampleName=$(basename $bam | sed 's:.bam::g')

		for scaller in $(echo $SCALLER | tr "," "\n")
		do
			
			case $scaller in
			    gatk|GATK)
					echo "$(tput setaf 2)ccSNP: running gatk on $bam$(tput sgr0)"
					#GCYCLES is defined by ccSNP parameter
					gatkPipeline $REFERENCE $bam $PLOIDY $QUALITY $sampleName $GCYCLES
			    ;;
			    samtools|SAMTOOLS|SAMtools|bcftools|BCFtools)
					echo "$(tput setaf 2)ccSNP: running BCFtools on $bam$(tput sgr0)"

 					bcftools mpileup -f $REFERENCE -Ou -q 60 -Q $QUALITY $bam | 
 					bcftools call --ploidy $PLOIDY -mv -Ov > ${sampleName}.raw.vcf

 					bcftools view --include 'FMT/GT="1" && QUAL>=100 && DP>=10' ${sampleName}.raw.vcf  | 
 					bcftools view --types snps |
 					vt normalize -r $REFERENCE - | 
 					bcftools annotate --remove '^FORMAT/GT,^FORMAT/DP' > ${sampleName}.bcftools.vcf

 					#if [ "$NOIN" == "false" ];then #integration will occur
 					#	cleanName=$(bcftools query -l ${sampleName}.bcftools.vcf)
 					#	sed "s/${cleanName}/${sampleName}_bcftools/g" ${sampleName}.bcftools.vcf > tmp
 					#	rm ${sampleName}.bcftools.vcf && mv tmp ${sampleName}.bcftools.vcf
 					#fi

			    ;;
			    freebayes|FREEBAYES|fb|FB|freeb)
					echo "$(tput setaf 2)ccSNP: running Freebayes on $bam$(tput sgr0)"


					fasta_generate_regions.py ${REFERENCE}.fai $(awk '{print int($2/31)}' ${REFERENCE}.fai) > ref.txt
					SECONDS=0
					freebayes-parallel ref.txt $(nproc) -p $PLOIDY -C 10 --min-repeat-entropy 1 --strict-vcf -q $QUALITY -m 30 --min-coverage 10  -f $REFERENCE $bam > ${sampleName}.raw.vcf
					bcftools view --include 'FMT/GT="1" && QUAL>=100 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0' ${sampleName}.raw.vcf | 
					bcftools view --types snps | 
					vt normalize -r $REFERENCE - | 
					bcftools annotate --remove '^INFO/DP,^FORMAT/GT,^FORMAT/GL' > ${sampleName}.freebayes.vcf
					rm -f ref.txt

			    ;;
			esac
		done
		rm -f *.raw.vcf*
			
		if [ "$NOIN" == "false" ]; then

			VCFFILES=$(ls ${sampleName}*.vcf | tr '\n' ',' | sed ' s/,$//g')
			callcore $VCFFILES ${sampleName}.isec
		fi

	done

	cd ..
	export VCFFILES=$(ls $(pwd)/3-callvariant/*.vcf | tr '\n' ',' | sed ' s/,$//g')
	echo "$(tput setaf 2)ccSNP: call variants step done$(tput sgr0)"
}
