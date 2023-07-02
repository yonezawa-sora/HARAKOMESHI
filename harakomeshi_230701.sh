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

cat << EOS
RUN : success!
EOS
