function mergeReads {
set -e
    echo "$(tput setaf 2)ccSNP: Quality control reads step start$(tput sgr0)"

	READS=$1
	QUALITY=$2

	mkdir -p 0-qc
	cd 0-qc

		echo $READS | awk -F';' '{
				n=split($1,r1,",");
				split($2,r2,",");
				for(i=1;i<=n;i++){
					print r1[i]"\t"r2[i]
				}
			}' | while read R1 R2
		do
			sampleName=$(recognizeSampleName $(basename $R1) $(basename $R2))
			fastp -i $R1 -I $R2 --merge --include_unmerged --fix_mgi_id --merged_out ${sampleName}.fastq --trim_front1 0 --trim_tail1 0 --trim_front2 0 --trim_tail2 0 --correction --thread $(nproc)
			

			mv fastp.html ${sampleName}.html
		done
		rm -f fastp.json
		export READS=$(readlink -f $(ls -1 *.fastq) | tr '\n' ',' | sed "s/,$//g")
		PAIRED="false"
	cd ..
	echo "$(tput setaf 2)ccSNP: Quality control reads step done$(tput sgr0)"

}