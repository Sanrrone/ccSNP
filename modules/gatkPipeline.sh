#!/bin/bash
function gatkPipeline {


#usage: bash gatkPipeline REF.fa mybam.bam [--spark=true]
#make sure you have gatk (alias of java -jar gatk.x.x.x.x), bcftools, vt, samtools.

REF=$1
BAM=$2
PLOIDY=$3
QUALITY=$4
sampleName=$5
ITER=0

MemoryAvail=$(cat /proc/meminfo | grep "MemAvailable" | awk '{printf "%.0fG\n",$2/1000/1000}')

if [[ "$@" =~ "--spark=true" ]];then
	HC="HaplotypeCallerSpark --spark-master local[$(nproc)]"
	BR="BaseRecalibratorSpark --spark-master local[$(nproc)]"
	ABR="ApplyBQSRSpark --spark-master local[$(nproc)]"
else
	HC="HaplotypeCaller"
	BR="BaseRecalibrator"
	ABR="ApplyBQSR"
fi
re='^[0-9]+$'
if [[ $6 =~ $re ]]; then ITER=$6; fi
if [ "$QUALITY" == "" ];then QUALITY=30; fi
if [ "$PLOIDY" == "" ] || [ $((PLOIDY)) -lt 0 ] ;then PLOIDY=1; fi

GATKVERSION=$(command -v gatk | awk -F'/' '{gsub("gatk-","",$(NF-1));print $(NF-1)}') 
if [ "$(printf '%s\n' "4.1.4.0" "$GATKVERSION" | sort -V | head -n1)" == "4.1.4.0" ]; then CONTROVERSIALFLAG="-I"; else CONTROVERSIALFLAG="-F"; fi


 for i in $(seq 1 1 $ITER)
 do
	#step1 get first vairants vcf
	echo "gatk $HC -R $REF -I $BAM"
	gatk $HC -R $REF -I $BAM --base-quality-score-threshold $QUALITY --minimum-mapping-quality $QUALITY --read-filter AllowAllReadsReadFilter \
	--sample-ploidy $PLOIDY --min-pruning  3 --native-pair-hmm-threads $(nproc) -O snps.raw.vcf

	#step2 filter variants
	bcftools view --include 'FMT/GT="1" && QUAL>=100 && FMT/DP>=10' snps.raw.vcf  | 
	vt normalize -r $REF - | 
	bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' > ${sampleName}.vcf
	gatk IndexFeatureFile $CONTROVERSIALFLAG ${sampleName}.vcf

 	#step3 Quality Base Recalibration Table
 	
 	gatk $BR -I $BAM -R $REF --known-sites ${sampleName}.vcf -O recalibration.table

 	#step 4 recalibration
 	gatk $ABR -I $BAM -R $REF --bqsr-recal-file recalibration.table -O recal_cycle${i}_${sampleName}.bam
	
	BAM="recal_cycle${i}_${sampleName}.bam"
	#delete temporary files
	rm -f recal_cycle$((i-1))_${sampleName}.bam ${sampleName}.vcf ${sampleName}.vcf.idx snps.raw.vcf snps.raw.vcf.idx
 done


#BAM now have the corrected bam
gatk $HC -R $REF -I $BAM --base-quality-score-threshold $QUALITY --minimum-mapping-quality $QUALITY --read-filter AllowAllReadsReadFilter \
--sample-ploidy $PLOIDY --min-pruning  3 --native-pair-hmm-threads $(nproc) -O snps.raw.vcf

bcftools view --include 'FMT/GT="1" && QUAL>=100 && FMT/DP>=10' snps.raw.vcf  |  
bcftools view --types snps |
vt normalize -r $REF - |
bcftools annotate --remove '^INFO/TYPE,^INFO/DP,^FORMAT/GT,^FORMAT/DP,^FORMAT/GL' > ${sampleName}.gatk.vcf
rm -f *.raw.vcf* recal_cycle${i}_${sampleName}.* recalibration.table

}