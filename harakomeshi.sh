#!/bin/bash
set -xe
set -o pipefail

# メモ：このスクリプトは、米澤奏良(https://github.com/yonezawa-sora)が勉強用に書き直しているスクリプトです｡
# 大部分は http://dojineko.hateblo.jp/entry/2016/06/30/225113 から引用させていただきました。(ikraより)


####################

# オプション関連は以下から
PROGNAME="$( basename $0 )"

VERSION="v1.0"

# <1.Usage >
# function usage()は、ヘルプを表示する関数
function usage() {
  cat << EOS >&2   # 標準エラー出力にリダイレクト   
Harakomeshi ${VERSION} -RNAseq pipeline centered on Salmon for plants-
Usage: ${PROGNAME} experiment_table.csv species [options]
  args
    1.experiment matrix(csv)
    2.reference (rice) #現在はriceのみ対応
Options:
  -t, --threads NUM  Number of threads (default: 4)

EOS
  exit 1
}


# <2.version>
# function version()は、バージョンを表示する関数
function version() {
  cat << EOS >&2
ikra ${VERSION} -RNAseq pipeline centered on Salmon for plants-
EOS
  exit 1
}

# <3.Default value>
# デフォルト値を先に定義しておく(ikraより)
# ここで定義した変数は、オプションで上書きされる?
RUNINDOCKER=1
DOCKER=docker
THREADS=4

# <4.parse options>
# "PARAM=()"で変数の指定が可能

PARAM=()
for opt in "$@"; do # $@は引数の配列
    case "${opt}" in # case文でオプションを判別
        '-t'|'--threads' ) # -tまたは--threadsを指定した場合
            THREADS=4; shift # THREADSに4を代入して、shiftで$1を削除
            if [[ -n "$1" ]] && [[ ! "$1" =~ ^-+ ]]; then # $1が空でなく、かつ、$1がオプションでない場合
                THREADS="$1"; shift # THREADSに$1を代入して、$1を削除
            fi
            ;;
        '-o'|'--output' )
            if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
            fi
              OUTPUT_FILE="$2"
              shift 2
              ;;
        '-l'|'--log' )
            if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
            fi
              LOG_FILE="$2"
              shift 2
                ;;
            H_GEN_VER="$2"
            M_GEN_VER="$2"
            shift 2
            ;;

        '-h' | '--help' )
            usage
            ;;
        '-v' | '--version' )
            version
            ;;
        '-r' | '--remove' )
            IF_REMOVE_INTERMEDIATES=true ; shift
            ;;
        '--' | '-' )
            shift
            PARAM+=( "$@" )
            break
            ;;
        
        -* )
            echo "${PROGNAME}: illegal option -- '$( echo $1 | sed 's/^-*//' )'" 1>&2
            exit 1
            ;;
        * )
            if [[ -n "$1" ]] && [[ ! "$1" =~ ^-+ ]]; then # $1が空でなく、かつ、$1がオプションでない場合
                PARAM+=( "$1" ); shift # $1をPARAMに追加して、$1を削除
            fi
            ;;
    esac # case文の終了
done

# <5.Check parameters> 

# 必須のパラメータ(ikraでは2つ)が指定されているかチェックする
EX_MATRIX_FILE="${PARAM}"; PARAM=("${PARAM[@]:1}") # PARAMの先頭をEX_MATRIX_FILEに代入して、PARAMの先頭を削除
REF_SPECIES="${PARAM}"; PARAM=("${PARAM[@]:1}") # PARAMの先頭をREF_SPECIESに代入して、PARAMの先頭を削除

# -zとは､変数が空であるかどうかを判定する｡
[[ -z "${EX_MATRIX_FILE}" ]] && usage # -zでEX_MATRIX_FILEが空であるかどうかを判定し、空であればusageを表示
[[ -z "${REF_SPECIES}" ]] && usage

# 規定外のオプションがある場合にはusageを表示
if [[ -n "${PARAM[@]}" ]]; then
    usage
fi


# <6.make log file>
cat << EOS | tee -a ${LOG_FILE} # tee -aでログファイルに書き込み
harakomeshi ${VERSION} -RNAseq pipeline centered on Salmon for plants-
EOS

date >> ${LOG_FILE} # 現在時刻を表示
pwd >> ${LOG_FILE} # カレントディレクトリを表示
whoami >> ${LOG_FILE} # ユーザー名を表示
uname -n >> ${LOG_FILE} # ホスト名を表示


# <7.オプションテスト用>
# 結果を表示(オプションテスト用)
cat << EOS | column -t | tee -a ${LOG_FILE}
THREADS ${THREADS}
OUTPUT_FILE ${OUTPUT_FILE}
LOG_FILE ${LOG_FILE}
EOS

set -u # 未定義変数を使おうとしたらエラーを出力する

####################


# 実験テーブル.csv

# 十分大きなものにする。
MAXSIZE=25G
SRA_ROOT=$HOME/ncbi/public/sra

SCRIPT_DIR=$(cd $(dirname $0); pwd)

if [[ $REF_SPECIES = mouse ]]; then
  BASE_REF_TRANSCRIPT=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M${M_GEN_VER}
  REF_TRANSCRIPT=gencode.vM${M_GEN_VER}.transcripts.fa.gz
  if [ $IF_PC = false ]; then
    REF_TRANSCRIPT=gencode.vM${M_GEN_VER}.transcripts.fa.gz
  else
    REF_TRANSCRIPT=gencode.vM${M_GEN_VER}.pc_transcripts.fa.gz
  fi
  SALMON_INDEX=salmon_index_mouse
#   REF_GTF=gencode.vM${M_GEN_VER}.annotation.gtf.gz
  TX2SYMBOL=gencode.vM${M_GEN_VER}.metadata.MGI.gz

elif [[ $REF_SPECIES = human ]]; then
  BASE_REF_TRANSCRIPT=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${H_GEN_VER}
  # REF_TRANSCRIPT=gencode.v${H_GEN_VER}.pc_transcripts.fa.gz

  if [ $IF_PC = false ]; then
    REF_TRANSCRIPT=gencode.v${H_GEN_VER}.transcripts.fa.gz
  else
    REF_TRANSCRIPT=gencode.v${H_GEN_VER}.pc_transcripts.fa.gz
  fi

  SALMON_INDEX=salmon_index_human
#   REF_GTF=gencode.v${H_GEN_VER}.annotation.gtf.gz
  TX2SYMBOL=gencode.v${H_GEN_VER}.metadata.HGNC.gz

else
  echo No reference speice!
  exit
fi

# 変数の定義
COWSAY=cowsay
FASTQ_DUMP=fastq-dump
FASTERQ_DUMP=fasterq-dump
FASTP=fastp # fastpを追加
SALMON=salmon
RSCRIPT_TXIMPORT=Rscript
WGET=wget
PIGZ=pigz
TAR=tar

if [[ "$RUNINDOCKER" -eq "1" ]]; then
  echo "RUNNING IN DOCKER"
  # docker を走らせ終わったらコンテナを削除。(-rm)ホストディレクトリをコンテナにマウントする。(-v)

  if [[ $DOCKER = docker ]]; then
    DRUN="$DOCKER run  -u `id -u`:`id -g` --rm -v $PWD:/home -e HOME=/home --workdir /home "
  elif [[ $DOCKER = udocker ]]; then
    DRUN="$DOCKER run --rm -v $PWD:/home --workdir /home "
  fi

  SCRIPT_DIR=`dirname "$0"`
  #--user=biodocker

# docker image version check
# fastpも後で追加しようと考えています
  COWSAY_IMAGE=docker/whalesay
  SRA_TOOLKIT_IMAGE=quay.io/biocontainers/sra-tools:2.10.9--pl526haddd2b5_0
  SALMON_IMAGE=combinelab/salmon:1.4.0
  RSCRIPT_TXIMPORT_IMAGE=fjukstad/tximport
  WGET_IMAGE=fjukstad/tximport
  PIGZ_IMAGE=genevera/docker-pigz
  TAR_IMAGE=fjukstad/tximport


  $DOCKER pull $COWSAY_IMAGE
  $DOCKER pull $SRA_TOOLKIT_IMAGE
  $DOCKER pull $SALMON_IMAGE
  $DOCKER pull $RSCRIPT_TXIMPORT_IMAGE
  $DOCKER pull $PIGZ_IMAGE
  $DOCKER pull $TAR_IMAGE

# 変数の定義を上書き
  COWSAY="$DRUN $COWSAY_IMAGE $COWSAY"
  FASTQ_DUMP="$FASTQ_DUMP"
  FASTQC="$DRUN $FASTQC_IMAGE $FASTQC" 
  FASTQ_DUMP="$FASTQ_DUMP"
  FASTERQ_DUMP="$FASTERQ_DUMP"
  SALMON="$DRUN $SALMON_IMAGE $SALMON"
  RSCRIPT_TXIMPORT="$DRUN $RSCRIPT_TXIMPORT_IMAGE $RSCRIPT_TXIMPORT"
  WGET="$DRUN $WGET_IMAGE $WGET"
  PIGZ="$DRUN $PIGZ_IMAGE"
  TAR="$DRUN $TAR_IMAGE $TAR"

# docker run --rm -v $PWD:/data -v $PWD:/root/ncbi/public/sra --workdir /data -it inutano/sra-toolkit bash
else
  echo "RUNNING LOCAL"
fi

echo $EX_MATRIX_FILE
cat $EX_MATRIX_FILE

# tximport_R.Rを削除
if [[  -f "tximport_R.R" ]]; then # 既にある場合は削除
  rm tximport_R.R
fi

# # tximport_R.Rを取ってくる。
# cp $SCRIPT_DIR/tximport_R.R ./

# 2019/06/09 devv1.3 tximport_R.Rを埋め込み(ikraより)

cat << 'EOF' > tximport_R.R
#! /usr/bin/Rscript
library(tximport)
library(readr)
library(stringr)
# Rscript tximport_R.R gencode.vM19.metadata.MGI.gz Illumina_PE_SRR.csv output.tsv
args1 = commandArgs(trailingOnly=TRUE)[1]
args2 = commandArgs(trailingOnly=TRUE)[2]
args3 = commandArgs(trailingOnly=TRUE)[3]
tx2knownGene <- read_delim(args1, '\t', col_names = c('TXNAME', 'GENEID'))
exp.table <- read.csv(args2, row.names=NULL)
files.raw <- exp.table[,2]
# files.raw <- c("SE/test/ttt30.fq.gz", "SE/test/ttt2.fq.gz")
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


####################

# Process1:fasterq_dump
  # SE
  if [ $LAYOUT = SE ]; then
    # fasterq_dump
    if [[ ! -f "$SRR.fastq.gz" ]]; then
      $FASTERQ_DUMP $SRR --threads $THREADS --force -p
      # gzip $SRR.fastq
      $PIGZ $SRR.fastq
    else
      echo "ERROR: Check the sra-tools version"
    
    fi

    # fastqc
    if [[ ! -f "${SRR}_fastqc.zip" ]]; then
      $FASTQC -t $THREADS ${SRR}.fastq.gz
    fi

  # PE
  else
    # fasterq_dump
    if [[ ! -f "${SRR}_1.fastq.gz" ]]; then
      $FASTERQ_DUMP $SRR --split-files --threads $THREADS --force -p
      $PIGZ ${SRR}_1.fastq
      $PIGZ ${SRR}_2.fastq
    else
      echo "ERROR: Check the sra-tools version"
    fi

    # fastqc
    if [[ ! -f "${SRR}_1_fastqc.zip" ]]; then
      $FASTQC -t $THREADS ${SRR}${SUFFIX_PE_1}
      $FASTQC -t $THREADS ${SRR}${SUFFIX_PE_2}
    fi
  fi
done
fi

if [[ ! -f "multiqc_report_raw_reads.html" ]]; then
  $MULTIQC -n multiqc_report_raw_reads.html .
fi

# determin threads for trim galore.
# the sweet spot for TG is 4
if [ $THREADS -gt 4 ] ; then
  THREADS_TRIMGALORE=4
else
  THREADS_TRIMGALORE=$THREADS
fi


for i in `tail -n +2  $EX_MATRIX_FILE | tr -d '\r'` 
do
  if [ $IF_FASTQ = false ]; then 
    # fasterq_dump
    name=`echo $i | cut -d, -f1`
    SRR=`echo $i | cut -d, -f2`
    LAYOUT=`echo $i | cut -d, -f3`
    dirname_fq="./"
  else
    name=`echo $i | cut -d, -f1`
    fq=`echo $i | cut -d, -f2`
    LAYOUT=`echo $i | cut -d, -f3`
    fqname_ext="${fq##*/}"
    # echo $fqname_ext

    # ファイル名を取り出す（拡張子なし）
    # basename_fq="${fqname_ext%.*.*}"
    basename_fq=${fqname_ext}
    dirname_fq=`dirname $fq`
    dirname_fq=${dirname_fq}/
    SRR=${basename_fq}
  fi


  # trim_galore
  # SE
  if [ $LAYOUT = SE ]; then
    if [  -f "${dirname_fq}${SRR}.fq" ] && [ ! -f "${dirname_fq}${SRR}.fastq.gz" ]; then
      $PIGZ ${dirname_fq}${SRR}.fq
      ln -s ${SRR}.fq.gz ${dirname_fq}${SRR}.fastq.gz
    fi
    if [ -f "${dirname_fq}${SRR}.fastq" ] && [ ! -f "${dirname_fq}${SRR}.fastq.gz" ]; then
      $PIGZ ${dirname_fq}${SRR}.fastq
    fi
    if [ -f "${dirname_fq}${SRR}.fq.gz" ] && [ ! -f "${dirname_fq}${SRR}.fastq.gz" ]; then
      ln -s ${SRR}.fq.gz ${dirname_fq}${SRR}.fastq.gz
    fi

    if [[ ! -f "${dirname_fq}${SRR}_trimmed.fq.gz" ]]; then
      $TRIMGALORE --cores ${THREADS_TRIMGALORE} ${dirname_fq}${SRR}.fastq.gz
    fi

    # fastqc
    if [[ ! -f "${dirname_fq}${SRR}_trimmed_fastqc.zip" ]]; then
      $FASTQC -t $THREADS ${dirname_fq}${SRR}_trimmed.fq.gz
    fi

  # PE
  else
    if [ -f "${dirname_fq}${SRR}_1.fq" ] && [ ! -f "${dirname_fq}${SRR}_1.fastq.gz" ]; then
      ${PIGZ} ${dirname_fq}${SRR}_1.fq
      ${PIGZ} ${dirname_fq}${SRR}_2.fq
      ln -s ${SRR}_1.fq.gz ${dirname_fq}${SRR}_1.fastq.gz
      ln -s ${SRR}_2.fq.gz ${dirname_fq}${SRR}_2.fastq.gz
    fi
    if [ -f "${dirname_fq}${SRR}_1.fastq" ] && [ ! -f "${dirname_fq}${SRR}_1.fastq.gz" ]; then
      $PIGZ ${dirname_fq}${SRR}_1.fastq
      $PIGZ ${dirname_fq}${SRR}_2.fastq
    fi
    if [  -f "${dirname_fq}${SRR}.fq.gz" ] && [ ! -f "${dirname_fq}${SRR}_1.fastq.gz" ]; then
      cp ${dirname_fq}${SRR}.fq.gz ${dirname_fq}${SRR}.fastq.gz
    fi
    if [ -f "${dirname_fq}${SRR}${SUFFIX_PE_1}" ] && [ ! -f "${dirname_fq}${SRR}_1.fastq.gz" ]; then
      ln -s ${SRR}${SUFFIX_PE_1} ${dirname_fq}${SRR}_1.fastq.gz
      ln -s ${SRR}${SUFFIX_PE_2} ${dirname_fq}${SRR}_2.fastq.gz
    fi

    # trimmomatic
    if [[ ! -f "${dirname_fq}${SRR}_1_val_1.fq.gz" ]]; then
      $TRIMGALORE --cores ${THREADS_TRIMGALORE} \
      --paired ${dirname_fq}${SRR}_1.fastq.gz ${dirname_fq}${SRR}_2.fastq.gz \
      --output_dir ${dirname_fq}
    fi

    # fastqc
    if [[ ! -f "${dirname_fq}${SRR}_1_val_1_fastqc.zip" ]]; then
      $FASTQC -t $THREADS ${dirname_fq}${SRR}_1_val_1.fq.gz
      $FASTQC -t $THREADS ${dirname_fq}${SRR}_2_val_2.fq.gz
    fi
  fi
done

# download $REF_TRANSCRIPT
if [[ ! -f "$REF_TRANSCRIPT" ]]; then
  $WGET $BASE_REF_TRANSCRIPT/$REF_TRANSCRIPT
fi

################################

# instance salmon index
if [[ ! -d "$SALMON_INDEX" ]]; then
  $SALMON index --threads $THREADS --transcripts $REF_TRANSCRIPT --index $SALMON_INDEX -k 31 --gencode
fi

for i in `tail -n +2  $EX_MATRIX_FILE | tr -d '\r'`
do
  if [ $IF_FASTQ = false ]; then
    # fasterq_dump
    name=`echo $i | cut -d, -f1`
    SRR=`echo $i | cut -d, -f2`
    LAYOUT=`echo $i | cut -d, -f3`
    dirname_fq=""
  else
    name=`echo $i | cut -d, -f1`
    fq=`echo $i | cut -d, -f2`
    LAYOUT=`echo $i | cut -d, -f3`
    fqname_ext="${fq##*/}"
    # echo $fqname_ext

    # ファイル名を取り出す（拡張子なし）
    # basename_fq="${fqname_ext%.*.*}"
    basename_fq=${fqname_ext}
    dirname_fq=`dirname $fq`
    dirname_fq=${dirname_fq}/
    SRR=${basename_fq}
  fi

  # SE
  if [ $LAYOUT = SE ]; then
    if [[ ! -f "salmon_output_${SRR}/quant.sf" ]]; then
      mkdir salmon_output_${SRR}
      # libtype auto detection mode
      $SALMON quant -i $SALMON_INDEX \
      -l A \
      -r ${dirname_fq}${SRR}_trimmed.fq.gz \
      -p $THREADS \
      -o salmon_output_${SRR} \
      --gcBias \
      --validateMappings
  #       -g $REF_GTF
    fi

   # PE
  else
    if [[ ! -f "salmon_output_${SRR}/quant.sf" ]]; then
      mkdir salmon_output_${SRR}
      # libtype auto detection mode
      $SALMON quant -i $SALMON_INDEX \
      -l A \
      -1 ${dirname_fq}${SRR}_1_val_1.fq.gz \
      -2 ${dirname_fq}${SRR}_2_val_2.fq.gz \
      -p $THREADS \
      -o salmon_output_${SRR} \
      --gcBias \
      --validateMappings
  #       -g $REF_GTF
    fi
  fi
done

# multiqc
if [[ ! -f "multiqc_report.html" ]]; then
  $MULTIQC -n multiqc_report.html .
fi

# download $TX2SYMBOL
if [[ ! -f "$TX2SYMBOL" ]]; then
  $WGET $BASE_REF_TRANSCRIPT/$TX2SYMBOL
fi

# tximport
if [[ ! -f "$OUTPUT_FILE" ]]; then
  $RSCRIPT_TXIMPORT tximport_R.R $TX2SYMBOL $EX_MATRIX_FILE $OUTPUT_FILE
fi

# tximport
if [[  -f "tximport_R.R" ]]; then
  rm tximport_R.R
fi

# if [[ "$RUNINDOCKER" -eq "1" ]]; then
#
#   chmod 755 .
#
# fi

cat << EOS | tee -a ${LOG_FILE}
RUN : success!
EOS
