function bamprep {
set -e

echo "$(tput setaf 2)ccSNP: bam preparation step start $(tput sgr0)"

BAMFILES=$1
REFERENCE=$2
mkdir -p 2-bamprep
cd 2-bamprep
mem="600M"

for bam in $(echo "$BAMFILES" | tr ',' ' ') # we know ',' is the delimiter between bamfiles
do
   
   sampleName=$(basename $bam | sed 's:.bam::g') 
   echo "$(tput setaf 2)-- extracting only mapped reads on $bam $(tput sgr0)"
   
   samtools  view -b -@ $(nproc) -f 260 $bam > unmapped.bam
   #picard AddOrReplaceReadGroups I=mapped.bam O=fixed.bam RGID=$sampleName RGLB=foo RGPU=bar RGPL=illumina RGSM=$sampleName CREATE_INDEX=False
   echo "$(tput setaf 2)-- sorting, filtering clips, duplicated and fix mate reads on $bam$(tput sgr0)"
   
   samtools sort -n -l 0 -T /tmp --threads $(nproc) -m $mem mapped.bam | 
   samtools view -h | 
   samclip --max 10 --invert --ref $REFERENCE > ${sampleName}.bam


   samtools index ${sampleName}.bam
   #rm -f fixed.bam m*.bam
   rm -f mapped.bam

done

cd ..

export BAMFILES=$(ls -1 $(pwd)/2-bamprep/*.bam | tr '\n' ',' | sed "s/,$//g")

echo "$(tput setaf 2)ccSNP: bam preparation step done$(tput sgr0)"

}

