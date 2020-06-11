function qc {

    echo "$(tput setaf 2)ccSNP: Quality control reads step start$(tput sgr0)"

	local READS=$1
	local QUALITY=$2

	mkdir -p 0-qc
	cd 0-qc

		echo $READS | awk -F';' '{n=split($1,r1,",");split($2,r2,",");for(i=1;i<=n;i++){print r1[i]"\t"r2[i]}}' | while read R1 R2
		do
			sampleName=$(recognizeSampleName $(basename $R1) $(basename $R2))
			R1="$(readlink -f ../../$R1)"
			R2="$(readlink -f ../../$R2)"
			if [ "$R2" == "" ];then
				fastp -i $R1 --trim_front1 5 --trim_tail1 5 --cut_mean_quality $QUALITY --cut_right --average_qual $QUALITY --length_required 50 --thread $(nproc) --out1 ${sampleName}.qc
			else
				fastp -i $R1 -I $R2 --merge --include_unmerged --fix_mgi_id --merged_out ${sampleName}.qc --trim_front1 5 --trim_tail1 5 --trim_front2 5 --cut_right --trim_tail2 5 \
				--cut_mean_quality $QUALITY --average_qual $QUALITY --length_required 50 --correction --thread $(nproc) 

			fi
			mv fastp.html ${sampleName}.html
		done
		rm -f fastp.json
	cd ..
	echo "$(tput setaf 2)ccSNP: Quality control reads step done$(tput sgr0)"

}