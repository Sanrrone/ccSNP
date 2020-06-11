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



	for bam in $(echo "$BAMFILES" | tr ',' ' ')
	do
		sampleName=$(basename $bam | sed 's:.bam::g')

		for scaller in $(echo $SCALLER | tr "," "\n")
		do
			
			case $scaller in
			    gatk|GATK)
					echo "$(tput setaf 2)ccSNP: running gatk $bam$(tput sgr0)"
					#GCYCLES is defined by ccSNP parameter
					gatkPipeline $REFERENCE $bam $PLOIDY $QUALITY $sampleName $GCYCLES
			    ;;
			    samtools|SAMTOOLS|SAMtools|bcftools|BCFtools)
					echo "$(tput setaf 2)ccSNP: running BCFtools $bam$(tput sgr0)"

 					bcftools mpileup -f $REFERENCE -Ou -q 60 -Q $QUALITY $bam | 
 					bcftools call --ploidy $PLOIDY -mv -Ov > ${sampleName}.raw.vcf

 					bcftools view --include 'FMT/GT="1" && QUAL>=100 && DP>=10' ${sampleName}.raw.vcf  | 
 					bcftools view --types snps |
 					vt normalize -r $REFERENCE - | 
 					bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^FORMAT/GT,^FORMAT/DP,^FORMAT/GL' > ${sampleName}.bcftools.vcf

 					#if [ "$NOIN" == "false" ];then #integration will occur
 					#	cleanName=$(bcftools query -l ${sampleName}.bcftools.vcf)
 					#	sed "s/${cleanName}/${sampleName}_bcftools/g" ${sampleName}.bcftools.vcf > tmp
 					#	rm ${sampleName}.bcftools.vcf && mv tmp ${sampleName}.bcftools.vcf
 					#fi

			    ;;
			    freebayes|FREEBAYES|fb|FB|freeb)
					echo "$(tput setaf 2)ccSNP: running Freebayes $bam$(tput sgr0)"


					fasta_generate_regions.py ${REFERENCE}.fai $(awk '{print int($2/31)}' ${REFERENCE}.fai) > ref.txt
					SECONDS=0
					freebayes-parallel ref.txt $(nproc) -p $PLOIDY -C 10 --min-repeat-entropy 1 --strict-vcf -q $QUALITY -m 30 --min-coverage 10  -f $REFERENCE $bam > ${sampleName}.raw.vcf
					bcftools view --include 'FMT/GT="1" && QUAL>=100 && FMT/DP>=10 && (FMT/AO)/(FMT/DP)>=0' ${sampleName}.raw.vcf | 
					bcftools view --types snps | 
					vt normalize -r $REFERENCE - | 
					bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^FORMAT/GT,^FORMAT/DP,^FORMAT/GL' > ${sampleName}.freebayes.vcf
					rm -f ref.txt

			    ;;
			esac
		done
		rm -f *.raw.vcf*
			
		if [ "$NOIN" == "false" ]; then

			VCFFILES=$(ls ${sampleName}*.vcf | tr '\n' ',' | sed ' s/,$//g')
			callcore $VCFFILES ${sampleName}.isec
			#for vcf in $(echo $VCFFILES | tr ',' '\n')
			#do
			#	rm -f $vcf
			#done
			#ISEC=$(ls ${sampleName}.*.vcf | head -n 1)
			#sampleCallName=$(basename $ISEC | sed 's/.vcf//g')
			#bcftools view $ISEC -Ob -o ${sampleCallName}.bcf
			#bcftools index -f --threads "$(nproc)" ${sampleCallName}.bcf
			#rm $ISEC
			#ISEC="${sampleCallName}.bcf"
#
			#for vcf in $(ls ${sampleName}.*.vcf)
			#do
			#    sampleCallName=$(basename $vcf | sed 's/.vcf//g')
			#    bcftools view $vcf -Ob -o ${sampleCallName}.bcf
			#    bcftools index -f --threads "$(nproc)" ${sampleCallName}.bcf
			#    bcftools isec -Ob ${sampleCallName}.bcf $ISEC -p .
			#    mv 0002.bcf isec.bcf
			#    ISEC="isec.bcf"
			#    rm $vcf
			#done
			#rm -f ${sampleName}.*.bcf ${sampleName}.*.bcf.csi 000[0123].bcf 000[0123].bcf.csi
			#bcftools view $ISEC -Ov -o ${sampleName}.isec.vcf
			#rm -f $ISEC README.txt
		fi

	done

	cd ..
	export VCFFILES=$(ls $(pwd)/${OUTPUT}_3-callvariant/*.vcf | tr '\n' ',' | sed ' s/,$//g')
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
