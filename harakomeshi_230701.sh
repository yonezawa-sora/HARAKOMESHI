#!/bin/bash
set -xe
set -o pipefail

# derived from ikra 

PROGNAME="$( basename $0 )" # Program name

# Usage
function usage() {
    echo "Usage: $0 <SRR csv file> [threads [VALUE]]"
    echo "Options:"
    echo "  -h, --help    Show usage."
    exit 1
}

# Get Arguments
csv_file=$1

shift # shift command moves the positional parameters to the left by one!!

# Parse options （Referring to ikra.sh）
while (( $# )); do
    case "$1" in
        '-t'|'--threads' )
            THREADS=4; shift
            if [[ -n "$1" ]] && [[ ! "$1" =~ ^-+ ]]; then
                THREADS="$1"; shift
            fi
            ;;
        -h|--help)
            usage
            ;;
        *)
            break
            ;;
    esac
done

# Argument Checking
# '$# -ne 2' means 'if the number of arguments is not 2'

if [ -z "$THREADS" ]; then
    THREADS=4
fi

# Read from text file one row at a time
while IFS=, read -r name srr layout condition; do
    echo "Processing $name, $srr, $layout, $condition"

    # Separate process for SE and PE
    if [ "$layout" = "SE" ]; then
        # if SE...
        fastp -i "${srr}.fastq.gz" -o "${srr}_trimmed.fq.gz" -h "${srr}.html" -j "${srr}.json" -w $THREADS
        # If fastp command is successful, remove the original fastq file
        rm "${srr}.fastq.gz"

    elif [ "$layout" = "PE" ]; then
        # if PE...
        fastp -i "${srr}_1.fastq.gz" -I "${srr}_2.fastq.gz" -o "${srr}_1_trimmed.fq.gz" -O "${srr}_2_trimmed.fq.gz" -h "${srr}.html" -j "${srr}.json" -w $THREADS --detect_adapter_for_pe
        # If fastp command is successful, remove the original fastq files
        rm "${srr}_1.fastq.gz"
        rm "${srr}_2.fastq.gz"
    else
        echo "Invalid layout: $layout"
        echo "See fastp --help"
        exit 1
    fi 

done < "$csv_file"

cat << EOS
RUN : success!
EOS
