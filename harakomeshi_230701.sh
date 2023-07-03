#!/bin/bash
set -xe
set -o pipefail

# derived from ikra 

PROGNAME="$( basename $0 )" # Program name

VERSION="v0.0.1" # Version

# Usage
function usage() {
    cat << EOS >&2 
${PROGNAME} ${VERSION} 
    Usage: $0 <SRR csv file> [-t | --threads [VALUE], -h | --help]
    args:
        SRR csv file (ikra format)
    Optional args:
        -t, --threads: number of threads (default: 4)
        -h, --help: print help
EOS
    exit 1
}

# Default value

DRUN=`docker run -u $(id -u):$(id -g) --rm -v $(pwd):/home -e HOME=/home --workdir /home`
SALMON_INDEX='salmon_index_rice'
SALMON_IMAGE='combinelab/salmon:1.10.0'
GET_REF_TRANSCRIPTS='curl -O https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/oryza_sativa/cdna/Oryza_sativa.IRGSP-1.0.cdna.all.fa.gz'
REF_TRANSCRIPT='Oryza_sativa.IRGSP-1.0.cdna.all.fa.gz'
RSCRIPT_TXIMPORT=Rscript
RSCRIPT_TXIMPORT_IMAGE=fjukstad/tximport


# Get Arguments
csv_file=$1; shift # shift command moves the positional parameters to the left by one!!

# Parse options （Referring to ikra.sh）
PARAM=()
for opt in "$@"; do
    case "$opt" in
        '-t'|'--threads' )
            THREADS=4; shift
            if [[ -n "$1" ]] && [[ ! "$1" =~ ^-+ ]]; then
                THREADS="$1"; shift
            fi
            ;;
        '-h'| '--help' )
            usage
            ;;
        *)
            break
            ;;
    esac
done



# Arguments Check
if [ -z "$csv_file" ]; then # 引数が空の場合､usageを表示
    usage
fi

################################## fastp ##################################

tail -n +2 $csv_file | tr -d '\r' | while read i; do
    name=$(echo $i | cut -d ',' -f 1)
    SRR=$(echo $i | cut -d ',' -f 2)
    LAYOUT=$(echo $i | cut -d ',' -f 3)
    dirname_fq="./"

# Separate process for SE and PE
    if [ "$LAYOUT" = "SE" ]; then
    # if SE...
    fastp -i "${SRR}.fastq.gz" -o "${SRR}_trimmed.fq.gz" -h "${SRR}.html" -j "${SRR}.json" -w $THREADS
    # fastp command is successful, remove the original fastq file
    rm "${SRR}.fastq.gz"

    elif [ "$LAYOUT" = "PE" ]; then
     # if PE...
    fastp -i "${SRR}_1.fastq.gz" -I "${SRR}_2.fastq.gz" -o "${SRR}_1_trimmed.fq.gz" -O "${SRR}_2_trimmed.fq.gz" -h "${SRR}.html" -j "${SRR}.json" -w $THREADS --detect_adapter_for_pe
    # fastp command is successful, remove the original fastq files
    rm "${SRR}_1.fastq.gz"
    rm "${SRR}_2.fastq.gz"
    else
    echo "Invalid layout: $LAYOUT"
    echo "See fastp --help"
    exit 1
    fi 
done

################################## Salmon ##################################

# docker run
SALMON="$DRUN $SALMON_IMAGE $SALMON"

# instance salmon index

if [[ ! -d $SALMON_INDEX ]]; then
    $GET_REF_TRANSCRIPTS
    eval $SALMON index \ 
    --threads $THREADS --transcripts $REF_TRANSCRIPT --index $SALMON_INDEX -k 31
fi

# quantification

tail -n +2 $csv_file | tr -d '\r' | while read i; do
    name=$(echo $i | cut -d ',' -f 1)
    SRR=$(echo $i | cut -d ',' -f 2)
    LAYOUT=$(echo $i | cut -d ',' -f 3)
    dirname_fq=""

# SE
    if [ $LAYOUT = "SE" ]; then
        if [[ ! -f "salmon_output_${SRR}/quant.sf" ]]; then
        mkdir salmon_output_${SRR}
        # libtype auto detection mode (ikra) 
        $SALMON quant -i $SALMON_INDEX \
        -l A \
        -r ./${SRR}_trimmed.fq.gz \
        -p $THREADS \
        -o salmon_output_${SRR} \
        --gcbias \
        --validateMappings
        fi
    # PE
    else
        if [[ ! -f "salmon_output_${SRR}/quant.sf" ]]; then
        mkdir salmon_output_${SRR}
        # libtype auto detection mode (ikra)
        $SALMON quant -i $SALMON_INDEX \
        -l A \
        -1 ./${SRR}_1_trimmed.fq.gz \
        -2 ./${SRR}_2_trimmed.fq.gz \
        -p $THREADS \
        -o salmon_output_${SRR} \
        --gcbias \
        --validateMappings
        fi
    fi
done

################################## tximport ##################################

# docker run
RSCRIPT_TXIMPORT="$DRUN $RSCRIPT_TXIMPORT_IMAGE $RSCRIPT_TXIMPORT"

# tximport script embedded in this script
cat << 'EOF' > tximport_R.R
#! /usr/bin/Rscript
library(tximport)
library(readr)
library(stringr)
args1 = commandArgs(trailingOnly=TRUE)[1]
args2 = commandArgs(trailingOnly=TRUE)[2]
args3 = commandArgs(trailingOnly=TRUE)[3]
tx2knownGene <- read_delim(args1, '\t', col_names = c('TXNAME', 'GENEID'))
exp.table <- read.csv(args2, row.names=NULL)
files.raw <- exp.table[,2]
files.raw <- gsub(".gz$", "", files.raw)
files.raw <- gsub(".fastq$", "", files.raw)
files.raw <- gsub(".fq$", "", files.raw)
split.vec <- sapply(files.raw, basename)
# print(paste(c("salmon_output_") , split.vec, c("/quant.sf"), sep=''))
# files <- paste(c("salmon_output_") , exp.table[,2], c("/quant.sf"), sep='')
files <- paste(c("salmon_output_") , split.vec, c("/quant.sf"), sep='')
names(files) <- exp.table[,1]
print(files)
# txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2knownGene)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2knownGene, countsFromAbundance="scaledTPM")
write.table(txi.salmon$counts, file=args3, sep="\t",col.names=NA,row.names=T,quote=F,append=F)
write.table(exp.table[-c(2,3)], file="designtable.csv",row.names=F,quote=F,append=F)
EOF

cat << EOS
RUN : success!
EOS
