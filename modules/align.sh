function align {
set -e

ALIGNER=$1
READS=$2
REFERENCE=$3
REFIDX=$(basename $REFERENCE | awk '{n=split($1,a,".");for(i=1;i<n;i++){if(i<n-1){printf("%s.",a[i])}else{printf("%s",a[i])}};printf "\n"}')

	echo "$(tput setaf 2)ccSNP: align step start$(tput sgr0)"
	mkdir -p 1-map
	################################### Mapping zone
	echo $READS | awk -F';' '{n=split($1,r1,",");split($2,r2,",");for(i=1;i<=n;i++){print r1[i]"\t"r2[i]}}' | while read R1 R2
	do
		if [ "$R2" == "" ];then
			sampleName=$(recognizeSampleName $(basename $R1))
		else
			sampleName=$(recognizeSampleName $(basename $R1) $(basename $R2))
		fi
		reads="$R1 $R2"
			cd 1-map
			    case $ALIGNER in
			        bwa|BWA)
						echo "$(tput setaf 2)BWA: running bwa mem on $reads$(tput sgr0)"
						bwa mem  -Y -M -R "@RG\tID:${sampleName}\tSM:${sampleName}\tPL:ILLUMINA\tPG:bwa" -t $(nproc) $REFERENCE $reads | samtools view -b -o ${sampleName}.bwa.bam -@ $(nproc)
			        ;;
			        smalt|SMALT)
 						smalt map -f sam -o ${sampleName}.smalt.bam -q 10 -n $(nproc) $REFERENCE $reads #minimum quality set to 10
			        ;;
			        novoalign|NOVOALIGN)
						novoalign -f $REFERENCE -d $ref -c $(nproc) -o SAM | samtools view -b -o ${sampleName}.novoalign.bam -@ $(nproc)
			        ;;
			    esac
		    cd ..
	done



export BAMFILES=$(ls -1 $(pwd)/1-map/*.bam | tr '\n' ',')
echo "$(tput setaf 2)ccSNP: align step done$(tput sgr0)"

}