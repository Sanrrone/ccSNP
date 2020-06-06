function callcore {

VCFFILES=$1

for vcf in $(echo $VCFFILES)
do
    sampleName=$(basename $vcf | sed 's/.vcf//g')
    bcftools view $vcf -Ob -o ${sampleName}.bcf
    bcftools index --threads "$(nproc)" ${sampleName}.bcf

done
bcftools merge -Ov -o merged.vcf *.bcf
awk -F'\t' '{for(i=9;i<=NF;i++){if($i==".:."){next}};print}' merged.vcf > coreSNP.vcf
rm -f *.bcf *.csi merged.vcf


}