#!/bin/bash 

# ヒト、マウスの解析の部分とtximportの部分を除いてみました（220905）

set -xe

# basename コマンドを使ってスクリプト名だけ取り出している
# $0は起動中のシェルスクリプト名がセットされている

PROGNAME="$( basename $0 )"

# []の外は必須のものとして認識される []はオプションを意味する
# speciesは文字列として認識させる
# Usage
function usage() {
  cat << EOS >&2        
 -RNAseq pipeline centered on Salmon-
Usage: ${PROGNAME} experiment_table.csv [--test --fastq, --help, --without-docker, --udocker, --protein-coding] [--transcriptome [VALUE]][--genome [VALUE]][--threads [VALUE]][--output [VALUE]][--suffix_PE_1 [VALUE]][--suffix_PE_2 [VALUE]]
  args
    1.experiment matrix(csv)
    2.reference(human or mouse)
Options:
  --test  test mode(MAX_SPOT_ID=100000).(dafault : False)
  --fastq use fastq files instead of SRRid. The extension must be foo.fastq.gz (default : False)
  -u, --udocker
  -w, --without-docker
  -pc, --protein-coding use protein coding transcripts instead of comprehensive transcripts. (defalut : True)
  -ct, --comprehensive-transcripts use comprehensive transcripts instead of protein coding transcripts. (default : False) 
  -t, --threads
  -o, --output  output file. (default : output.tsv)  
  -l, --log  log file. (default : ikra.log)
  -a, --align carry out mapping onto a reference genome. hisat2 or star (default : None)
  -g, --gencode specify the version of gencode. (defalut : Mouse=26, Human=37)
  -s1, --suffix_PE_1    suffix for PE fastq files. (default : _1.fastq.gz)
  -s2, --suffix_PE_2    suffix for PE fastq files. (default : _2.fastq.gz)
  -h, --help    Show usage.
  -v, --version Show version.
  -r, --remove-intermediates Remove intermediate files


Github repo : https://github.com/yyoshiaki/ikra
EOS
  exit 1
}

# デフォルト値を先に定義しておく
RUNINDOCKER=1
DOCKER=docker
THREADS=1
IF_TEST=false
IF_FASTQ=false
IF_PC=True
SUFFIX_PE_1=_1.fastq.gz
SUFFIX_PE_2=_2.fastq.gz
OUTPUT_FILE=output.tsv
LOG_FILE=ikra.log
MAPPING_TOOL=None
IF_REMOVE_INTERMEDIATES=false

#オプションをパース
#"$@は引数一つ一つがそれぞれ別のものとして認識される"
#"PARAM=()"で変数の指定が可能


PARAM=()
for opt in "$@"; do
    case "${opt}" in
        #　モード選択など引数の無いオプションの場合
        '--test' )
            IF_TEST=true; shift
            ;;
        '--fastq' )
            IF_FASTQ=true; shift
            ;;
        '-pc'|'--protein-coding' )
            IF_PC=true; shift
            ;;
        '-ct'|'--comprehensive-transcripts' )
            IF_PC=true; shift
            ;;
        '-u'|'--udocker' )
            DOCKER=udocker; shift
            ;;
        '-w'|'--without-docker' )
            RUNINDOCKER=0; shift
            ;;
        #　引数が任意の場合
        '--transcriptome' )
        # outputfileの部分をコピーして変数を変更して適用
            if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
            fi
              REF_TRANSCRIPT="$2"
              shift 2
              ;; 

        '--genome' )
           if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
            fi
              REF_GENOME="$2"
              shift 2
              ;;  

        '-t'|'--threads' )
            THREADS=4; shift
            if [[ -n "$1" ]] && [[ ! "$1" =~ ^-+ ]]; then
                THREADS="$1"; shift
            fi
            ;;
        '-s1'|'--suffix_PE_1' )
            if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
            fi
            SUFFIX_PE_1="$2"
            shift 2
            ;;

        '-s2'|'--suffix_PE_2' )
            if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
            fi
            SUFFIX_PE_2="$2"
            shift 2
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

        '-a'|'--align' )
            if [[ "$2" == "hisat2" ]]; then
                MAPPING_TOOL=HISAT2
            elif [[ "$2" == "star" ]]; then
                MAPPING_TOOL=STAR
            elif [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
            fi
              shift 2
                ;;

        '-g'|'--gencode' )
            if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
                exit 1
            fi
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
            if [[ -n "$1" ]] && [[ ! "$1" =~ ^-+ ]]; then
                PARAM+=( "$1" ); shift
            fi
            ;;
    esac
done


# オプション無しの値を使う場合はここで処理する
EX_MATRIX_FILE="${PARAM}"; PARAM=("${PARAM[@]:1}")
REF_SPECIES="${PARAM}"; PARAM=("${PARAM[@]:1}")

[[ -z "${EX_MATRIX_FILE}" ]] && usage
[[ -z "${REF_SPECIES}" ]] && usage

# 規定外のオプションがある場合にはusageを表示
if [[ -n "${PARAM[@]}" ]]; then
    usage
fi


# オプション無しの値を使う場合はここで処理する
EX_MATRIX_FILE="${PARAM}"; PARAM=("${PARAM[@]:1}")
REF_SPECIES="${PARAM}"; PARAM=("${PARAM[@]:1}")

[[ -z "${EX_MATRIX_FILE}" ]] && usage
[[ -z "${REF_SPECIES}" ]] && usage

# 規定外のオプションがある場合にはusageを表示
if [[ -n "${PARAM[@]}" ]]; then
    usage
fi


date >> ${LOG_FILE}
pwd >> ${LOG_FILE}
whoami >> ${LOG_FILE}
uname -n >> ${LOG_FILE}


# 結果を表示(オプションテスト用)
cat << EOS | column -t | tee -a ${LOG_FILE}
EX_MATRIX_FILE ${EX_MATRIX_FILE}
REF_SPECIES ${REF_SPECIES}
RUNINDOCKER ${RUNINDOCKER}
DOCKER ${DOCKER}
THREADS ${THREADS}
IF_TEST ${IF_TEST:-false}
IF_FASTQ ${IF_FASTQ:-false}
IF_PC ${IF_PC:-false}
IF_REMOVE_INTERMEDIATES ${IF_REMOVE_INTERMEDIATES:-false}
OUTPUT_FILE ${OUTPUT_FILE}
MAPPING_TOOL ${MAPPING_TOOL}
M_GEN_VER ${M_GEN_VER}
H_GEN_VER ${H_GEN_VER}
LOG_FILE ${LOG_FILE}
EOS

set -u



# オプション関連ここまで
# 十分大きなサイズにする
MAXSIZE=25G
SRA_ROOT=$HOME/ncbi/public/sra

# dirnameでシェルスクリプトを抜き出す
SCRIPT_DIR=$(cd $(dirname $0); pwd)


# リファレンス
SALMON_INDEX=salmon_index


# 変数の指定
COWSAY=cowsay
# PREFETCH=prefetch
FASTQ_DUMP=fastq-dump
FASTERQ_DUMP=fasterq-dump
FASTQC=fastqc
MULTIQC=multiqc
# TRIMMOMATIC=trimmomatic
TRIMGALORE=trim_galore
HISAT2=hisat2
STAR_MAPPING=STAR
SAMBAMBA=sambamba
BAMCOVERAGE=bamCoverage
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

  # 危険！
  # chmod 777 .


  COWSAY_IMAGE=docker/whalesay
  # quay.io/biocontainers/sra-tools:2.10.7--pl526haddd2b5_1 had an error.
  # the earlier version may stop during the download.
  SRA_TOOLKIT_IMAGE=quay.io/biocontainers/sra-tools:2.10.9--pl526haddd2b5_0
  FASTQC_IMAGE=biocontainers/fastqc:v0.11.9_cv8
  MULTIQC_IMAGE=quay.io/biocontainers/multiqc:1.10.1--py_0
#   TRIMMOMATIC_IMAGE=fjukstad/trimmomatic
#   TRIMMOMATIC_IMAGR=comics/trimmomatic
  TRIMGALORE_IMAGE=quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0
  HISAT2_IMAGE=quay.io/biocontainers/hisat2:2.2.1--h1b792b2_3
  STAR_IMAGE=quay.io/biocontainers/star:2.7.8a--h9ee0642_1
  SAMBAMBA_IMAGE=quay.io/biocontainers/sambamba:0.8.0--h984e79f_0
#   salmonのイメージ取得を1.4.0から1.9.0に変更  
  SALMON_IMAGE=combinelab/salmon:1.9.0
#   SALMON_IMAGE=fjukstad/salmon
  RSCRIPT_TXIMPORT_IMAGE=fjukstad/tximport
  WGET_IMAGE=fjukstad/tximport
  PIGZ_IMAGE=genevera/docker-pigz
  TAR_IMAGE=fjukstad/tximport
  BAMCOVERAGE_IMAGE=quay.io/biocontainers/deeptools:3.5.1--py_0


# docker pullする
  $DOCKER pull $COWSAY_IMAGE
  $DOCKER pull $SRA_TOOLKIT_IMAGE
  $DOCKER pull $FASTQC_IMAGE
  $DOCKER pull $MULTIQC_IMAGE
  # $DOCKER pull $TRIMMOMATIC_IMAGE
  $DOCKER pull $TRIMGALORE_IMAGE
  $DOCKER pull $HISAT2_IMAGE
  $DOCKER pull $STAR_IMAGE
  $DOCKER pull $SAMBAMBA_IMAGE
  $DOCKER pull $BAMCOVERAGE_IMAGE
  $DOCKER pull $SALMON_IMAGE
  $DOCKER pull $RSCRIPT_TXIMPORT_IMAGE
  $DOCKER pull $PIGZ_IMAGE
  $DOCKER pull $TAR_IMAGE


  COWSAY="$DRUN $COWSAY_IMAGE $COWSAY"
  # PREFETCH="$DRUN -v $PWD:/root/ncbi/public/sra $SRA_TOOLKIT_IMAGE $PREFETCH"
  # FASTQ_DUMP="$DRUN $SRA_TOOLKIT_IMAGE $FASTQ_DUMP"
  FASTQ_DUMP="$FASTQ_DUMP"
#  FASTERQ_DUMP="$DRUN $SRA_TOOLKIT_IMAGE $FASTERQ_DUMP"
  FASTQC="$DRUN $FASTQC_IMAGE $FASTQC" 
  FASTQ_DUMP="$FASTQ_DUMP"
  FASTERQ_DUMP="$FASTERQ_DUMP"
  MULTIQC="$DRUN $MULTIQC_IMAGE $MULTIQC"
#   TRIMMOMATIC="$DRUN $TRIMMOMATIC_IMAGE $TRIMMOMATIC"
  # TRIMMOMATIC="$DRUN $TRIMMOMATIC_IMAGE " # fjukstad/trimmomaticのentrypointのため
  TRIMGALORE="$DRUN $TRIMGALORE_IMAGE $TRIMGALORE"
  HISAT2="$DRUN $HISAT2_IMAGE $HISAT2"
  STAR_MAPPING="$DRUN $STAR_IMAGE $STAR_MAPPING"
  SAMBAMBA="$DRUN $SAMBAMBA_IMAGE $SAMBAMBA"
  BAMCOVERAGE="$DRUN $BAMCOVERAGE_IMAGE $BAMCOVERAGE"
  SALMON="$DRUN $SALMON_IMAGE $SALMON"
#   SALMON="$DRUN $SALMON_IMAGE"
  RSCRIPT_TXIMPORT="$DRUN $RSCRIPT_TXIMPORT_IMAGE $RSCRIPT_TXIMPORT"
  WGET="$DRUN $WGET_IMAGE $WGET"
  PIGZ="$DRUN $PIGZ_IMAGE"
  TAR="$DRUN $TAR_IMAGE $TAR"


   # docker run --rm -v $PWD:/data -v $PWD:/root/ncbi/public/sra --workdir /data -it inutano/sra-toolkit bash
else
  echo "RUNNING LOCAL"
fi



# if [ $MAX_SPOT_ID = 0 ]; then
if [ $IF_TEST = true ]; then
  $COWSAY "test mode( MAX_SPOT_ID is set)"
  MAX_SPOT_ID="-X 100000"
else
  MAX_SPOT_ID=""
fi


echo $EX_MATRIX_FILE
cat $EX_MATRIX_FILE

# tximport
# 現段階ではsalmonのindex作って定量できるかまで確認する(220906)



# fasterq_dump
  # SE
  if [ $LAYOUT = SE ]; then
    # fastq_dump
    if [[ ! -f "$SRR.fastq.gz" ]]; then
      if [[ $MAX_SPOT_ID == "" ]]; then
        $FASTERQ_DUMP $SRR --threads $THREADS --force -p
        # gzip $SRR.fastq
        $PIGZ $SRR.fastq
      else
        $FASTQ_DUMP $SRR $MAX_SPOT_ID --gzip
      fi
    fi

    # fastqc
    if [[ ! -f "${SRR}_fastqc.zip" ]]; then
      $FASTQC -t $THREADS ${SRR}.fastq.gz
    fi

  # PE
  else
    # fastq_dump
    if [[ ! -f "${SRR}_1.fastq.gz" ]]; then
      if [[ $MAX_SPOT_ID == "" ]]; then
        $FASTERQ_DUMP $SRR --split-files --threads $THREADS --force -p
        # gzip ${SRR}_1.fastq
        # gzip ${SRR}_2.fastq
        $PIGZ ${SRR}_1.fastq
        $PIGZ ${SRR}_2.fastq
      else
        $FASTQ_DUMP $SRR $MAX_SPOT_ID --gzip --split-files
      fi
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


# download $REF_TRANSCRIPTは今回行わない
# mapping by hisat2

################################


# instance salmon index
if [[ ! -d "$SALMON_INDEX" ]]; then
  $SALMON index --threads $THREADS --transcripts $REF_TRANSCRIPT --index $SALMON_INDEX -k 31 
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


if [ $IF_REMOVE_INTERMEDIATES = true ]; then
  rm -f *fastq.gz
  rm -f *fq.gz
  rm -f *fastqc.zip
  rm -rf salmon_output_*
fi
# if [[ "$RUNINDOCKER" -eq "1" ]]; then
#
#   chmod 755 .
#
# fi

cat << EOS | tee -a ${LOG_FILE}
RUN : success!
EOS