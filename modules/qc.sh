function qc {
set -e
    echo "$(tput setaf 2)ccSNP: Quality control reads step start$(tput sgr0)"

	READS=$1
	QUALITY=$2

	mkdir -p 0-qc
	cd 0-qc

		echo $READS | awk -F';' '{
			if($2==""){
				n=split($1,r0,",");
				for(i=1;i<=n;i++){
					print r0[i]
				}
			}else{
				n=split($1,r1,",");
				split($2,r2,",");
				for(i=1;i<=n;i++){
					print r1[i]"\t"r2[i]
				}
			}}' | while read R1 R2
		do
			if [ "$R2" == "" ];then
				sampleName=$(recognizeSampleName $(basename $R1))
				fastp -i $R1 --trim_front1 5 --trim_tail1 5 --cut_mean_quality $QUALITY --cut_right --average_qual $QUALITY --length_required 30 --thread $(nproc) --out1 ${sampleName}.qc.fastq
			else
				sampleName=$(recognizeSampleName $(basename $R1) $(basename $R2))
				fastp -i $R1 -I $R2 --merge --include_unmerged --fix_mgi_id --merged_out ${sampleName}.qc.fastq --trim_front1 5 --trim_tail1 5 --trim_front2 5 --cut_right --trim_tail2 5 \
				--cut_mean_quality $QUALITY --average_qual $QUALITY --length_required 40 --correction --thread $(nproc) 
			fi

			mv fastp.html ${sampleName}.html
		done
		rm -f fastp.json
		export READS=$(readlink -f $(ls -1 *.qc.fastq) | tr '\n' ',' | sed "s/,$//g")
		PAIRED="false"
	cd ..
	echo "$(tput setaf 2)ccSNP: Quality control reads step done$(tput sgr0)"

}