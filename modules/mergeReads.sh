function mergeReads {
set -e
    echo "$(tput setaf 2)ccSNP: Quality control reads step start$(tput sgr0)"

	READS=$1

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
			fastp -i $R1 -I $R2 --merge --include_unmerged --fix_mgi_id --merged_out ${sampleName}.fastq --fix_mgi_id --correction --thread $(nproc) -Q 
			mv fastp.html ${sampleName}.html
		done
		rm -f fastp.json
		export READS=$(readlink -f $(ls -1 *.fastq) | tr '\n' ',' | sed "s/,$//g")
		PAIRED="false"
	cd ..
	echo "$(tput setaf 2)ccSNP: Quality control reads step done$(tput sgr0)"

}