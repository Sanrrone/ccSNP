function bamprep {
set -e

echo "$(tput setaf 2)ccSNP: bam preparation step start $(tput sgr0)"

BAMFILES=$1
REFERENCE=$2
OUTPUT=$3
mkdir -p ${OUTPUT}_2-bamprep
cd ${OUTPUT}_2-bamprep
mem="600M"

for bam in $(echo "$BAMFILES" | tr ',' ' ') # we know ',' is the delimiter between bamfiles
do
   
   sampleName=$(basename $bam | sed 's:.bam::g') 
   #mapped within the insert size and in correct orientation
   #samtools view -b -@ $(nproc) -f 99 $bam > m1.bam
   #samtools view -b -@ $(nproc) -f 147 $bam > m2.bam
   #samtools view -b -@ $(nproc) -f 83 $bam > m3.bam
   #samtools view -b -@ $(nproc) -f 163 $bam > m4.bam

   #if [ PAIRED="true" ];then
   #   #One of the reads is unmapped
   #   samtools view -b -@ $(nproc) -f 73 $bam > m5.bam
   #   samtools view -b -@ $(nproc) -f 133 $bam > m6.bam
   #   samtools view -b -@ $(nproc) -f 89 $bam > m7.bam
   #   samtools view -b -@ $(nproc) -f 121 $bam > m8.bam
   #   samtools view -b -@ $(nproc) -f 165 $bam > m9.bam
   #   samtools view -b -@ $(nproc) -f 181 $bam > m10.bam
   #   samtools view -b -@ $(nproc) -f 101 $bam > m11.bam
   #   samtools view -b -@ $(nproc) -f 117 $bam > m12.bam
   #   samtools view -b -@ $(nproc) -f 153 $bam > m13.bam
   #   samtools view -b -@ $(nproc) -f 185 $bam > m14.bam
   #   samtools view -b -@ $(nproc) -f 69 $bam > m15.bam
   #   samtools view -b -@ $(nproc) -f 137 $bam > m16.bam
   #fi
   #echo "merging bams"
   #samtools merge mapped.bam m*.bam
   echo "$(tput setaf 2)ccSNP: extracting only mapped reads on $bam $(tput sgr0)"
   
   samtools  view -b -@ $(nproc) -F 260 $bam > mapped.bam
   #picard AddOrReplaceReadGroups I=mapped.bam O=fixed.bam RGID=$sampleName RGLB=foo RGPU=bar RGPL=illumina RGSM=$sampleName CREATE_INDEX=False
   echo "$(tput setaf 2)ccSNP: sorting, filtering clips, duplicated and fix mate reads$(tput sgr0)"
   
   samtools sort -n -l 0 -T /tmp --threads $(nproc) -m $mem mapped.bam | 
   samtools view -h | 
   samclip --max 10 --ref $REFERENCE | 
   samtools fixmate -m - - | 
   samtools sort -l 0 -T /tmp --threads $(nproc) -m $mem | 
   samtools markdup -T /tmp -r -s - - > ${sampleName}.filt.bam


   samtools index ${sampleName}.filt.bam
   #rm -f fixed.bam m*.bam
   rm -f mapped.bam
done

cd ..

export BAMFILES=$(ls -1 $(pwd)/${OUTPUT}_2-bamprep/*.filt.bam | tr '\n' ',' | sed ' s/,$//g')

echo "$(tput setaf 2)ccSNP: bam preparation step done$(tput sgr0)"

}

