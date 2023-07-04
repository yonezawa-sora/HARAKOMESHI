#!/bin/bash
set -xe

# derived from ikra 

PROGNAME="$( basename $0 )" # Program name

VERSION="v1.0.0" # Version

# Usage
function usage() {
    cat << EOS >&2 
${PROGNAME} ${VERSION} 
    Usage: $0 <SRR csv file> <output.tsv file> [-t | --threads [VALUE], -h | --help]
    args:
        SRR csv file (ikra format)
        output.tsv file (for tximport)
    Optional args:
        -t, --threads: number of threads (default: 4)
        -h, --help: print help
EOS
    exit 1
}

# Default value
# バッククォート(``)ではコマンド置換を意味するため､変数の定義ではなく､実際に実行してしまう
# そのため､ ダブルクォート("")で囲むことで､変数の定義として扱う

DRUN="docker run -u $(id -u):$(id -g) --rm -v $(pwd):/home -e HOME=/home --workdir /home"
SALMON_INDEX="salmon_index_rice"
SALMON_IMAGE="combinelab/salmon:1.10.0"
GET_REF_TRANSCRIPTS="curl -O https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/oryza_sativa/cdna/Oryza_sativa.IRGSP-1.0.cdna.all.fa.gz"
REF_TRANSCRIPT="Oryza_sativa.IRGSP-1.0.cdna.all.fa.gz"
RSCRIPT_TXIMPORT_IMAGE="fjukstad/tximport"
TX2GENEID="rice_tx2geneID.tsv"


# Get Arguments
csv_file=$1; shift # shift command moves the positional parameters to the left by one!!
output_file=$1; shift

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

# 最後の文字が改行文字(\n)でなくても､読み込めるように変更
tail -n +2 $csv_file | tr -d '\r' | while IFS= read -r i || [[ -n "$i" ]]; do
    name=$(echo $i | cut -d ',' -f 1)
    SRR=$(echo $i | cut -d ',' -f 2)
    LAYOUT=$(echo $i | cut -d ',' -f 3)
    dirname_fq="./"

# Separate process for SE and PE
    if [ "$LAYOUT" = "SE" ]; then
    # if SE...
        if [[ ! -f "${SRR}_trimmed.fq.gz" ]]; then
        fastp -i "${SRR}.fastq.gz" -o "${SRR}_trimmed.fq.gz" -h "${SRR}.html" -j "${SRR}.json" -w $THREADS
        # fastp command is successful, remove the original fastq file
        rm "${SRR}.fastq.gz"
        fi

    elif [ "$LAYOUT" = "PE" ]; then
     # if PE...
        if [[ ! -f "${SRR}_1_trimmed.fq.gz" ]]; then
        fastp -i "${SRR}_1.fastq.gz" -I "${SRR}_2.fastq.gz" -o "${SRR}_1_trimmed.fq.gz" -O "${SRR}_2_trimmed.fq.gz" -h "${SRR}.html" -j "${SRR}.json" -w $THREADS --detect_adapter_for_pe
        # fastp command is successful, remove the original fastq files
        rm "${SRR}_1.fastq.gz"
        rm "${SRR}_2.fastq.gz"
        fi
    else
    echo "Invalid layout: $LAYOUT"
    echo "See fastp --help"
    exit 1
    fi 
done

################################## Salmon ##################################

# instance salmon index

if [[ ! -d $SALMON_INDEX ]]; then
    if [[ ! -f $REF_TRANSCRIPT ]]; then
        $GET_REF_TRANSCRIPTS
    fi
    $DRUN $SALMON_IMAGE salmon index \
    --threads $THREADS --transcripts $REF_TRANSCRIPT --index $SALMON_INDEX -k 31
fi

# quantification

tail -n +2 $csv_file | tr -d '\r' | while IFS= read -r i || [[ -n "$i" ]]; do
    name=$(echo $i | cut -d ',' -f 1)
    SRR=$(echo $i | cut -d ',' -f 2)
    LAYOUT=$(echo $i | cut -d ',' -f 3)
    dirname_fq=""

# SE 
# remove option --gcBias for salmon newest version
    if [ $LAYOUT = "SE" ]; then
        if [[ ! -f "salmon_output_${SRR}/quant.sf" ]]; then
        mkdir salmon_output_${SRR}
        # libtype auto detection mode (ikra) 
        $DRUN $SALMON_IMAGE salmon quant -i $SALMON_INDEX \
        -l A \
        -r ./${SRR}_trimmed.fq.gz \
        -p $THREADS \
        -o salmon_output_${SRR} \
        --gcBias \
        --validateMappings
        fi
    # PE
    # remove option --gcBias for salmon newest version (1.10.0)
    else
        if [[ ! -f "salmon_output_${SRR}/quant.sf" ]]; then
        mkdir salmon_output_${SRR}
        # libtype auto detection mode (ikra)
        $DRUN $SALMON_IMAGE salmon quant -i $SALMON_INDEX \
        -l A \
        -1 ./${SRR}_1_trimmed.fq.gz \
        -2 ./${SRR}_2_trimmed.fq.gz \
        -p $THREADS \
        -o salmon_output_${SRR} \
        --gcBias \
        --validateMappings
        fi
    fi

    if [[ -f "salmon_output_${SRR}/logs/salmon_quant.log" ]]; then
    mv "salmon_output_${SRR}/logs/salmon_quant.log" "salmon_output_${SRR}/logs/${SRR}_salmon_quant.log"
    fi

done

################################## tximport ##################################

# tximport script embedded in this script(ikra)
cat << 'EOF' > tximport_R.R
#! /usr/bin/Rscript
##### <import 3 libraries> #####
library(tximport)
library(readr)
library(stringr)

##### <get 3 arguments([1]=TX2GENEID, [2]=$CSV_FILE, [3]=$OUTPUT_FILE)> #####
args1 = commandArgs(trailingOnly=TRUE)[1]
args2 = commandArgs(trailingOnly=TRUE)[2]
args3 = commandArgs(trailingOnly=TRUE)[3]

##### <Storage TX2GENEID as data.frame> #####
tx2knownGene <- read_delim(args1, '\t', col_names = c('TXNAME', 'GENEID'))

##### <Read CSV_FILE into a data frame> #####
exp.table <- read.csv(args2, row.names=NULL)

##### <Remove specific extensions from the file names in the second column> #####
files.raw <- exp.table[,2]
files.raw <- gsub(".gz$", "", files.raw)
files.raw <- gsub(".fastq$", "", files.raw)
files.raw <- gsub(".fq$", "", files.raw)

##### <Extract only the file names from the paths> #####
split.vec <- sapply(files.raw, basename)

##### <Strage file vector as list> #####
files <- paste(c("salmon_output_") , split.vec, c("/quant.sf"), sep='')
names(files) <- exp.table[,1]

##### <Print file path vector> #####
print(files)

##### <tximport process> #####
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2knownGene, countsFromAbundance="scaledTPM")

##### <Write txi.salmon$counts> #####
write.table(txi.salmon$counts, file=args3, sep="\t",col.names=NA,row.names=T,quote=F,append=F)
write.table(exp.table[-c(2,3)], file="designtable.csv",row.names=F,quote=F,append=F)
EOF

# tximport

if [[ ! -f "$output_file" ]]; then
    $DRUN $RSCRIPT_TXIMPORT_IMAGE Rscript ./tximport_R.R $TX2GENEID $csv_file $output_file
fi

if [[  -f "tximport_R.R" ]]; then
  rm tximport_R.R
fi

cat << EOS
RUN : success!
EOS
