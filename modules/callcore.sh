function callcore {

VCFFILES=$1
OUTPUTNAME=$2

if [ "$OUTPUTNAME" == "" ];then
	OUTPUTNAME="coreSNP"
fi 

for vcf in $(echo $VCFFILES | tr ',' '\n')
do
    sampleName=$(basename $vcf | sed 's/.vcf//g')
    bcftools view $vcf -Ob -o ${sampleName}.bcf
    bcftools index -f --threads "$(nproc)" ${sampleName}.bcf
done
bcftools merge -Ov -o merged.vcf *.bcf
awk -F'\t' '{for(i=9;i<=NF;i++){if($i==".:."){next}};print}' merged.vcf > ${OUTPUTNAME}.vcf
rm -f *.bcf *.csi merged.vcf

}

function consensus {

VCFFILES=$1
REFERENCE=$2

	for vcf in $(ls *.vcf)
	do
		sampleName=$(basename $vcf | sed 's/.vcf//g')
		bcftools view $vcf -Oz -o ${sampleName}.bcf
		bcftools index -f --threads "$(nproc)" ${sampleName}.bcf
		REFERENCE=$(echo $sampleName | awk -F"_" '{gsub("ref","",$1);print $1".fasta"}')
		echo "consensus on $REFERENCE - $sampleName"
		bcftools consensus -f $REFERENCE --sample ${sampleName} ${sampleName}.bcf > consensus_${sampleName}.fna
		sed -i "s/>.*/>${sampleName}/" consensus_${sampleName}.fna
	done

}