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

if [ $(echo "$READS" | tr ',' '\n' | wc -l | awk '{print $1}' ) -le 1 ];then
    NOCC="true"
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
if [ "$FORCE" == "" ];then
    FORCE=false
fi

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

###### check for WGET
if ! [ -x "$(command -v wget)" ]; then
    curl -O https://ftp.gnu.org/gnu/wget/wget2-1.99.2.tar.gz
    tar xzf wget2-1.99.2.tar.gz
    cd wget2-1.99.2
        ./configure --prefix=$EXTBINARIES/wget2-1.99
        make -j $(nproc)
        make install
    cd ..
    rm -rf wget2-1.99.2.tar.gz wget2-1.99.2
    alias wget=$EXTBINARIES/wget2-1.99/bin/wget2
    alias wget >/dev/null 2>&1 && 
    echo "$(tput setaf 2)wget found in $(command -v wget)$(tput sgr0)" \
    || \
    echo "$(tput setaf 1)ERROR: wget not installed$(tput sgr0)" && exit 1
else
    echo "$(tput setaf 2)wget found in $(command -v wget)$(tput sgr0)"
fi

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
    SCALLER="gatk" #default
    NOIN="true" #default
else 
    if [ "$SCALLER" == "all" ];then
        SCALLER="freebayes,samtools,gatk" #default
        NOIN="false" #default
    fi
fi


if [ $(echo "$SCALLER" | tr ',' '\n' | wc -l) -le 1 ];then
    NOIN="true"
else
    NOIN="false"
fi

for scaller in $(echo $SCALLER | tr "," "\n")
do
    case $scaller in
        gatk|GATK)
            if ! [ -x "$(command -v gatk)" ]; then
                wget https://github.com/broadinstitute/gatk/releases/download/4.1.7.0/gatk-4.1.7.0.zip
                unzip gatk-4.1.7.0.zip
                mv gatk-4.1.7.0 $EXTBINARIES/.
                rm gatk-4.1.7.0.zip
                alias gatk=$EXTBINARIES/gatk-4.1.7.0/gatk


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
        samtools|SAMTOOLS|SAMtools|bcftools|BCFtools)
            if [ -x "$(command -v bcftools)" ]; then
                 echo "$(tput setaf 2)bcftools found in $(command -v bcftools)$(tput sgr0)"
            fi
        ;;
        snver|SNVer|SNVER)
            if ! [ -x "$(command -v snver)" ]; then

                wget -O snver.tar.gz https://sourceforge.net/projects/snver/files/latest/download
                mkdir -p snver
                tar xzf snver.tar.gz -C snver
                rm -f snver.tar.gz
                mv snver $EXTBINARIES/.
                alias snver="java -jar $EXTBINARIES/snver/SNVerIndividual.jar"
                
                alias snver >/dev/null 2>&1 && 
                    echo "$(tput setaf 2)snver found in $(command -v snver)$(tput sgr0)" \
                    || \
                    echo "$(tput setaf 1)ERROR: snver not installed$(tput sgr0)" && exit 1
            else
                 echo "$(tput setaf 2)snver found in $(command -v snver)$(tput sgr0)"
            fi
        ;;
        freebayes|FREEBAYES|fb|FB|freeb)
            if ! [ -x "$(command -v freebayes)" ]; then

                if ! [ -x "$(command -v cmake)" ]; then
                    echo "$(tput setaf 1)cmake not found, you can avoid this installation installing directly freebayes$(tput sgr0)"
                    exit 1
                fi

                echo "$(tput setaf 2)ccSNP: downloading and installing freebayes and their all dependencies (I wish it works)$(tput sgr0)"

                wget https://github.com/ekg/freebayes/releases/download/v1.3.1/freebayes-v1.3.1
                chmod +x freebayes-v1.3.1
                mv freebayes-v1.3.1 $EXTBINARIES/freebayes

                wget https://raw.githubusercontent.com/ekg/freebayes/master/scripts/fasta_generate_regions.py
                chmod +x fasta_generate_regions.py
                mv fasta_generate_regions.py $EXTBINARIES/.
                
                wget https://raw.githubusercontent.com/ekg/freebayes/master/scripts/freebayes-parallel
                chmod +x freebayes-parallel
                mv freebayes-parallel $EXTBINARIES/.

                wget https://ftp.gnu.org/gnu/parallel/parallel-20200522.tar.bz2
                tar xjf parallel-20200522.tar.bz2
                cd parallel-20200522
                    ./configure --prefix=$EXTBINARIES/parallel-20200522
                    make -j $(nproc)
                    make install
                cd ..
                rm -rf parallel-20200522*
                cp $EXTBINARIES/parallel-20200522/bin/parallel $EXTBINARIES/.
                #alias parallel=$EXTBINARIES/parallel-20200522/bin/parallel


                git clone --recursive https://github.com/vcflib/vcflib.git
                cd vcflib
                    make -j $(nproc)
                    cp bin/* $EXTBINARIES/.
                cd ..
                rm -rf vcflib
                wget https://raw.githubusercontent.com/vcflib/vcflib/master/scripts/vcffirstheader
                chmod +x vcffirstheader
                mv vcffirstheader $EXTBINARIES/.


                if ! [ -x "$(command -v freebayes)" ]; then
                    echo "$(tput setaf 1)ERROR: freebayes is not installed." >&2
                    exit 1
                else
                    echo "$(tput setaf 2)freebayes found in $(command -v freebayes)$(tput sgr0)"
                fi
            else
                if ! [ -x "$(command -v fasta_generate_regions.py)" ]; then
                    wget https://raw.githubusercontent.com/ekg/freebayes/master/scripts/fasta_generate_regions.py
                    chmod +x fasta_generate_regions.py
                    mv fasta_generate_regions.py $EXTBINARIES/.
                fi
                if ! [ -x "$(command -v freebayes-parallel)" ]; then
                    wget https://raw.githubusercontent.com/ekg/freebayes/master/scripts/freebayes-parallel
                    chmod +x freebayes-parallel
                    mv freebayes-parallel $EXTBINARIES/.

                    wget https://ftp.gnu.org/gnu/parallel/parallel-20200522.tar.bz2
                    tar xjf parallel-20200522.tar.bz2
                    cd parallel-20200522
                        ./configure --prefix=$EXTBINARIES/parallel-20200522
                        make -j $(nproc)
                        make install
                    cd ..
                    rm -rf parallel-20200522*
                    alias parallel=$EXTBINARIES/parallel-20200522/bin/parallel

                    git clone --recursive https://github.com/vcflib/vcflib.git
                    cd vcflib
                        make -j $(nproc)
                        cp bin/* $EXTBINARIES/.
                    cd ..
                    rm -rf vcflib
                else
                    echo "$(tput setaf 2)freebayes found in $(command -v freebayes)$(tput sgr0)"
                fi
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