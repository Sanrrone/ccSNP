function callcore {

VCFFILES=$1
OUTPUTNAME=$2

if [ "$OUTPUTNAME" == "" ];then
	OUTPUTNAME="coreSNP"
fi 

#for vcf in $(echo $VCFFILES | tr ',' '\n')
#do
#    sampleName=$(basename $vcf | sed 's/.vcf//g')
#    bcftools view $vcf -Ob -o ${sampleName}.bcf
#    bcftools index -f --threads "$(nproc)" ${sampleName}.bcf
#done
#bcftools merge -Ov -o merged.vcf *.bcf
#awk -F'\t' '{for(i=9;i<=NF;i++){if($i==".:."){next}};print}' merged.vcf > ${OUTPUTNAME}.vcf
##mv merged.vcf ${OUTPUTNAME}_all.vcf
#rm -f *.bcf *.csi merged.vcf

ISEC=$(echo $VCFFILES | tr ',' '\n' | head -n 1)
sampleName=$(basename $ISEC | sed 's/.vcf//g')
cleanName=$(bcftools query -l $ISEC)
bcftools view $ISEC -Ob -o ${sampleName}.bcf
bcftools index -f --threads "$(nproc)" ${sampleName}.bcf
if [ "$NOIN" == "false" ]; then rm $ISEC; fi


ISEC="${sampleName}.bcf"

for vcf in $(echo $VCFFILES | tr ',' '\n' |awk '{if(NR>1)print}')
do
    sampleName=$(basename $vcf | sed 's/.vcf//g')
    bcftools view $vcf -Ob -o ${sampleName}.bcf
    bcftools index -f --threads "$(nproc)" ${sampleName}.bcf
    bcftools isec -Ob $ISEC ${sampleName}.bcf -p .
    mv 0002.bcf isec.bcf
    mv 0002.bcf.csi isec.bcf.csi
    ISEC="isec.bcf"
    if [ "$NOIN" == "false" ]; then rm $vcf; fi
done
rm -f ${sampleName}.*.bcf ${sampleName}.*.bcf.csi 000[0123].bcf 000[0123].bcf.csi
bcftools view $ISEC -Ov -o tmp
rm -f $ISEC README.txt *.bcf *.bcf.csi

sed "s/${cleanName}/${OUTPUTNAME}/g" tmp > ${OUTPUTNAME}.vcf
rm -f tmp

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