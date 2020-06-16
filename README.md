# ccSNP
Call Core SNP is a pipeline written in `bash` to combine three different variant caller softwares: BCFtools, Freebayes, GATK4. The pipeline take your reads and a reference as input and give you a VCF file with all shared snps between the samples.

# Try the test before install something
ccSNP will trying to download all the necessary binaries from their sources so the only requisite you need is have installed:

* git
* Java
* cmake
* Curl
* libssl-dev
* Python >= 3.7 (setted as python binary)

Try run the example and check if it runs without errors. You should have the coreSNP file inside the ccsnp folder and the Elephant in ascii :D.

Other way try to install by yourselfe the requisites.

# Requisites

Make sure you have these programs in your PATH variable:

* Samtools >= v1.7
* BCFtools >= v1.7
* Freebayes >=v1.3.1
* GATK >= v4.1.1.3.0
* BWA >= v0.7
* fastp >= v0.20.1(optional if you want to do QC step to your reads)

# Usage

##### Simple paired end reads
`ccSNP -1 reads_R1.fastq -2 reads_R2.fastq -r reference.fasta`
##### Simple single end reads
`ccSNP -0 reads.fastq -r reference.fasta`
##### Simple paired end reads with quality score snp filter of 30 and choosing bcftools and freebayes as variant caller.
`ccSNP -1 reads_R1.fastq -2 reads_R2.fastq -r reference.fasta -q 30 -c bcftools,freebayes`
##### Simple paired end reads using GATK as unique variant caller and 2 GATK cycles (base quality recalibration, 0 by default)
`ccSNP -1 reads_R1.fastq -2 reads_R2.fastq -r reference.fasta -c gatk -gc 2`
##### Multiple samples paired end reads and perform QC filter on the reads
`ccSNP -1 sample1_R1.fastq,sample2_R1.fastq -2 sample1_R2.fastq,sample2_R2.fastq -r reference.fasta -mqc`
##### Multiple samples single end reads with all variant callers and keep tmp files (the bams used, the vcf from each variant call, and the intersections)
`ccSNP -0 sample1.fastq,sample2.fastq,sample3.fastq -r reference.fasta -c all --noclean`



## Available options:

* `-1/-2`: for paired end reads, multiple samples can be added separated with ','.
* `-0`: for single reads, multiple samples can be added separated with ','.
* `-r`: reference file in fasta format.
* `-c`: varian Caller to use. By default it will use only bcftools, the options are bcftools, freebayes and gatk. you can select all of some of them separating the option with comma.
* `-o`: Output folder.
* `-q`: Quality score for quality filtering steps (QC reads and SNP filters). By default: 20.
* `-qc`: Perform quality filtering and trimmed reads usinq by default Q>=20 (see -q).
* `-gc`: GATK cycles, The best practices for GATK actually include the base quality recalibration and then the SNP call, with `-gc` flag you can set how much cycles you want. By default is 0.
* `--noin`: No Integrate results. By default ccSNP will intersect the results of a sample with diffrent snp callaers, with this option all the results for different callers and same samples are keept and not intersected.
* `--nocc`: No call core. By default ccSNP will calculate the core SNP between all samples, with this option there is no core SNP.
* `--noclean`: By default ccSNP delete results not used in the final steps like bam files of map step, sample vcf files before the intersection, tmp folders. Set this flag to keep all files
* `--force`: by default if the output folder of the same previous run exist the pipeline stop unless you give the --force flag, this will overwrite files with the same name.
* `--debug`: show all logs from all steps.


# Notes

* You can use one, two or three variant callers, is not necessary install all but the intersection of three results is better than only one.
* ccSNP only deal with SNP. Indels and other type of variant are filtered.
* Ploidy is set to 1 for all cases
* Filter for post SNP call is performed using DP>=10 and QUAL>=100