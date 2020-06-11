# ccSNP
Call Core SNP is a pipeline to combine three different variant caller softwares: BCFtools, Freebayes, GATK4. The pipeline take your reads and a reference as input and give you a VCF file with all shared snps between the samples.

# Try the test before install something
ccSNP will trying to download all the necessary binaries from their sources so the only requisite you need is have installed:

* git
* Java
* cmake
* Curl

Try run the example and check if it runs without errors. You should have the coreSNP file inside the ccsnp folder.

Other way try to install by yourselfe the requisites.

# Requisites

Make sure you have these programs in your PATH variable:
* Samtools >= v1.7
* BCFtools >= v1.7
* Freebayes
* GATK >= v4.1.1.3.0
* BWA >= v0.7

# Usage

* `ccSNP -1 reads_R2.fastq -2 reads_R2.fastq -r reference.fasta`
* `ccSNP -0 reads.fastq -r reference.fasta`
* `ccSNP -1 reads_R2.fastq -2 reads_R2.fastq -r reference.fasta -q 20 -c bcftools,freebayes`
* `ccSNP -1 reads_R2.fastq -2 reads_R2.fastq -r reference.fasta -c gatk`

## Available options:

* `-1/-2`: for paired end reads, multiple samples can be added separated with ',' e.g. ccSNP -1 sample1_r1.fastq,sample2_r1.fastq -2 sample1_r2.fastq,sample2_r2.fastq
* `-0`: for single reads, multiple samples can be added separated with ',' e.g. ccSNP -0 sample1.fastq,sample2.fastq,sample3.fastq -r reference.fasta
* `-r`: reference file in fasta format.
* `-c`: varian Caller to use. By default it will use only bcftools, the options are bcftools, freebayes and gatk. you can select all of some of them separating the option with comma. e.g. ccSNP -0 sample.fastq -r references.fasta -c bcftools,freebayes.
* `-o`: Output folder.
* `--noin`: No Integrate results. By default ccSNP will intersect the results of a sample with diffrent snp callaers, with this option all the results for different callers and same samples are keept and not intersected.
* `--nocc`: No call core. By default ccSNP will calculate the core SNP between all samples, with this option there is no core SNP.
* `--force`: by default if the output folder of the same previous run exist the pipeline stop unless you give the --force flag, this will overwrite files with the same name.
* `--debug`: show all logs from all steps.

# Example

Simple Paired end reads sample

`./ccSNP -1 example/EC958_R1.fastq -2 example/EC958_R2.fastq -r example/APECO1.fasta`

----

Multiple paired end reads samples

`./ccSNP -1 example/EC958_R1.fastq,example/MS6573_R1.fastq -2 example/EC958_R2.fastq,example/MS6573_R2.fastq -r example/APECO1.fasta`

----

Multiple paired end reads samples with multiple variant call

`./ccSNP -1 example/EC958_R1.fastq,example/MS6573_R1.fastq -2 example/EC958_R2.fastq,example/MS6573_R2.fastq -r example/APECO1.fasta -c samtools,freebayes,gatk`

