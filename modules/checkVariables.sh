function checkVariables {
	set -e
    echo "$(tput setaf 2)ccSNP: checking variables and files$(tput sgr0)"
####################################################################################################################
#paired reads check
if [ "$PAIRED" == "" ]; then
    PAIRED="true"

    if [ "$R1" == "" ] || [ "$R2" == "" ];then
        echo "$(tput setaf 1)ERROR: forward or reverse flags are empty, specify the file using -1/--forward and -2/--reverse flags"
        exit 1
    fi
##### check if R1 file exist
    for r1 in $(echo $R1 | tr "," "\n")
    do
      if [ ! -f "$r1" ];then
        echo "$(tput setaf 1)ERROR: Read $r1 doesn't exist in $(pwd), aborting"
        exit 1
      fi
    done
##### check if R2 file exist
    for r2 in $(echo $R2 | tr "," "\n")
    do
      if [ ! -f "$r2" ];then
        echo "$(tput setaf 1)ERROR: Read $r2 doesn't exist in $(pwd), aborting"
        exit 1
      fi
    done
    export READS="$R1;$R2"

    if [ $(echo $READS | awk -F';' '{n1=split($1,r1,",");n2=split($2,r2,",");if(n1 == n2){print "true"}else{print "falsa"}}') == "false" ];then
        echo "$(tput setaf 1)ERROR: paired reads doesnt match"
        exit 1
    fi

fi
# single read check
if [ "$PAIRED" == "false" ]; then
    if [ "$R0" == "" ];then
        echo "$(tput setaf 1)ERROR: unpaired read flag specified but parameter is empty"
        exit 1
    fi
    for r0 in $(echo $R0 | tr "," "\n")
    do
      if [ ! -f "$r0" ];then
        echo "$(tput setaf 1)ERROR: Read $r0 doesn't exist, aborting"
        exit 1
      fi
    done

    export READS="$R0"
fi


###################################################################################################################
#check quality variable

if [ "$QUALITY" == "" ];then
    QUALITY=20
else
    if [ $((QUALITY)) lt 5 ];then
        echo "$(tput setaf 1)ERROR: quality value too low (< 10)"
        exit 1
    fi
fi

##check output folder
##########################
if [ "$OUTPUT" == "" ];then
    OUTPUT="ccsnp"
fi

if [ -d "$OUTPUT" ];then
    if [ "$FORCE" == "false" ];then
        echo "$(tput setaf 1)ERROR: $OUTPUT directory already exist or use --force to continue (overwriting)"
        exit 1   
    fi
else
    mkdir -p $OUTPUT
fi



###########################
#check ploidy
if [ "$PLOIDY" == "" ];then
    PLOIDY=1
else
    if [ $(echo "$1" | awk -F"\n" '{print ($0 != $0+0)}') ];then
        if [ $((PLOIDY)) -lt 1 ];then
            echo "$(tput setaf 1)ERROR: --ploidy must be a number > 0"
            exit 1
        fi
    else
        echo "$(tput setaf 1)ERROR: --ploidy must be a number > 0"
        exit 1
    fi
fi

#check reference
if [ "$REFERENCE" == "" ];then
    echo "$(tput setaf 1)ERROR: reference must be provided, use -r/--reference flag"
    exit 1
else
    if ! [ -f $REFERENCE ];then
        echo "$(tput setaf 1)ERROR: reference $REFERENCE not exist"
        exit 1
    else
        REFERENCE=$(readlink -f $REFERENCE)
    fi
fi


###################################################################################################################
####################check binaries

###check for aligners
if [ "$ALIGNER" == "" ];then
    ALIGNER="bwa" #default
    if ! [ -x "$(command -v bwa)" ]; then
        git clone https://github.com/lh3/bwa
        cd bwa && make -j 4
            chmod +x bwa
            cp bwa $EXTBINARIES/.
        cd ..
        rm -rf bwa

    fi
    if ! [ -x "$(command -v bwa)" ]; then
        echo '$(tput setaf 1)ERROR: bwa is not installed.' >&2
        exit 1
    fi
fi

for aligner in $(echo $ALIGNER | tr "," "\n")
do
    case $aligner in
        bwa|BWA)
            if ! [ -x "$(command -v bwa)" ]; then
                echo "$(tput setaf 1)ERROR: bwa is not installed." >&2
                exit 1
            else
                echo "$(tput setaf 2)bwa found in $(command -v bwa)$(tput sgr0)"
            fi
        ;;
        smalt|SMALT)
            if ! [ -x "$(command -v smalt)" ]; then
                echo "$(tput setaf 1)ERROR: smalt is not installed." >&2
                exit 1
            else
                echo "$(tput setaf 2)smalt found in $(command -v smalt)$(tput sgr0)"
            fi
        ;;
        novoalign|NOVOALIGN)
            if ! [ -x "$(command -v novoalign)" ]; then
                echo "$(tput setaf 1)ERROR: novoalign is not installed." >&2
                exit 1
            fi
        ;;
    esac
done

###check for SNP CALLERS
if [ "$SCALLER" == "" ];then
    SCALLER="gatk"
fi

for scaller in $(echo $SCALLER | tr "," "\n")
do
    case $scaller in
        gatk|GATK)
            if ! [ -x "$(command -v gatk)" ]; then
                wget https://github.com/broadinstitute/gatk/releases/download/4.1.7.0/gatk-4.1.7.0.zip
                unzip gatk-4.1.7.0.zip
                mv gatk-4.1.7.0 $EXTBINARIES/.
                alias gatk=$EXTBINARIES/gatk-4.1.7.0/gatk

                    if [ "$(whereis gatk)" == "" ]; then
                            echo "$(tput setaf 1)ERROR: gatk is not installed." >&2
                            exit 1
                        
                    fi
            else
                echo "$(tput setaf 2)gatk found in $(command -v gatk)$(tput sgr0)"
                if [ "GCYCLES" == "" ];then
                    GCYCLES=0
                fi
                if [ $((GCYCLES)) -lt 0 ];then
                    echo "$(tput setaf 1)ERROR: gatk cycles must be integer >= 0. e.g -gc 1" >&2
                    exit 1
                fi
                source $CCSNPHOME/modules/gatkPipeline.sh

            fi
        ;;
        samtools|SAMTOOLS|mpileup|SAMtools)
            if ! [ -x "$(command -v samtools)" ]; then
                echo "$(tput setaf 1)ERROR: samtools is not installed." >&2
                exit 1
            else
                 echo "$(tput setaf 2)samtools found in $(command -v samtools)$(tput sgr0)"
            fi
        ;;
        snver|SNVer|SNVER)
            if ! [ -x "$(command -v snver)" ]; then
                echo "$(tput setaf 1)ERROR: snver is not installed." >&2
                exit 1
            else
                 echo "$(tput setaf 2)snver found in $(command -v snver)$(tput sgr0)"
            fi
        ;;
        freebayes|FREEBAYES)
            if ! [ -x "$(command -v freebayes)" ]; then
                echo "$(tput setaf 1)ERROR: freebayes is not installed." >&2
                exit 1
            else
                 echo "$(tput setaf 2)freebayes found in $(command -v freebayes)$(tput sgr0)"
            fi
        ;;
    esac
done


#check fastp binary
if [ "$MAKEQC" == "" ];then
    MAKEQC=false
else

    if ! [ -x "$(command -v fastp)" ]; then
        wget http://opengene.org/fastp/fastp
        chmod a+x ./fastp
        mv fastp $EXTBINARIES/.

        if ! [ -x "$(command -v fastp)" ]; then
            echo "$(tput setaf 1)ERROR: fastp is not installed." >&2
            exit 1
        fi
    else
        echo "$(tput setaf 2)fastp found in $(command -v fastp) $(tput sgr0)"
    fi
    MAKEQC=true
fi

#check samclip
if ! [ -x "$(command -v samclip)" ]; then
    git clone https://github.com/tseemann/samclip.git
    cp samclip/samclip $EXTBINARIES/.
    rm -rf samclip

    if ! [ -x "$(command -v samclip)" ]; then
        echo "$(tput setaf 1)ERROR: samclip is not installed." >&2
        exit 1
    fi
else
     echo "$(tput setaf 2)samclip found in $(command -v samclip) $(tput sgr0)"
fi

#check samtools
if ! [ -x "$(command -v samtools)" ]; then
    git clone git://github.com/samtools/htslib.git
    git clone https://github.com/samtools/samtools.git
    cd samtools
        autoconf -Wno-syntax && ./configure
        make -j $(nproc)
        mv samtools $EXTBINARIES/.
    cd ..
    rm -rf htslib samtools

    if ! [ -x "$(command -v samtools)" ]; then
        echo "$(tput setaf 1)ERROR: samtools is not installed." >&2
        exit 1
    fi
else
     echo "$(tput setaf 2)samtools found in $(command -v samtools) $(tput sgr0)"
fi

#check bcftools
if ! [ -x "$(command -v bcftools)" ]; then
    git clone git://github.com/samtools/htslib.git
    git clone git://github.com/samtools/bcftools.git
    cd bcftools
        # The following is optional:
        autoheader && autoconf && ./configure --enable-libgsl --enable-perl-filters
        make -j $(nproc)
        mv bcftools $EXTBINARIES/.
    cd ..
    rm -rf htslib bcftools

    if ! [ -x "$(command -v bcftools)" ]; then
        echo "$(tput setaf 1)ERROR: bcftools is not installed." >&2
        exit 1
    fi
else
     echo "$(tput setaf 2)bcftools found in $(command -v bcftools) $(tput sgr0)"
fi

#check vt
if ! [ -x "$(command -v vt)" ];then
    wget https://github.com/atks/vt/archive/0.57721.tar.gz
    tar xzf 0.57721.tar.gz
    cd vt-0.57721
        make -j $(nproc)
        mv vt $EXTBINARIES/.
    cd ..
    rm -rf vt-0.57721 xzf 0.57721.tar.gz

    if ! [ -x "$(command -v vt)" ]; then
        echo "$(tput setaf 1)ERROR: vt is not installed. Install it manually from https://github.com/atks/vt" >&2
        exit 1
    fi
else
     echo "$(tput setaf 2)vt found in $(command -v vt)$(tput sgr0)"
fi

##check picard
#if ! [ -x "$(command -v picard)" ]; then
#    if [ -f $EXTBINARIES/picard.jar ];then 
#        echo "$(tput setaf 2)picard found in $EXTBINARIES/picard.jar$(tput sgr0)"
#        alias picard="java -jar $EXTBINARIES/picard.jar"
#    else
#        wget https://github.com/broadinstitute/picard/releases/download/2.22.9/picard.jar
#        mv picard.jar $EXTBINARIES/.
#        alias picard="java -jar $EXTBINARIES/picard.jar"
#    fi
#
#    if [ "$(whereis picard.jar)" == "" ]; then
#        echo "$(tput setaf 1)ERROR: picard is not installed." >&2
#        exit 1
#    fi
#fi

echo "$(tput setaf 2)ccSNP: all variables and files looks fine$(tput sgr0)"

}

function recognizeSampleName {
    #R1 ACT AS R1 AND R0 IF IS SINGLE READ
    R1=$1
    R2=$2

    if [ "$PAIRED" == "true" ];then
        sampleName=$(echo $R1";"$R2 | awk -F';' '{n1=split($1,r1,"");n2=split($2,r2,"");for(i=1;i<=n1;i++){if(r1[i]==r2[i]){printf("%s",r1[i])}else{break}}printf "\n"}')
        sampleName=$(echo $sampleName | sed "s/_[rR]$//g")
    else
        sampleName=$(echo $R1 | sed "s/.fastq//g" | sed "s/.gz//g" | sed "s/.zip//g" | sed "s/.bz2//g")
    fi

    echo "$sampleName"

}

function printHelp {
    echo "Usage: ccSNP -1 reads_R2.fastq -2 reads_R2.fastq -r reference.fasta"
    echo "Usage: ccSNP -0 reads.fastq -r reference.fasta"
    echo "Usage: ccSNP -1 reads_R2.fastq -2 reads_R2.fastq -r reference.fasta -q 20 -ploidy 1"
    echo -e "\nAvailable options:\n"

    exit
}