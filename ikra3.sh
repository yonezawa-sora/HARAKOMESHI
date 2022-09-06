#!/bin/bash
set -xe

# オプション関連ここから
# 大部分は http://dojineko.hateblo.jp/entry/2016/06/30/225113 から引用させていただきました。

# 変数 EX_MATRIX_FILE, REF_SPECIES はここで定義
# if [[ $IF_TEST = true ]]; then でテストモード用の実行が可能

# 今まで$1 = EX_MATRIX_FILEだったのを変更している
# 以降の$1をEX_MATRIX_FILEで置き換える必要がある？(必要なら修正お願いします...)


PROGNAME="$( basename $0 )"

VERSION="v2.0.1"

cat << "EOF" 
    __                       
 __/\ \                      
/\_\ \ \/'\   _ __    __     
\/\ \ \ , <  /\`'__\/'__`\   
 \ \ \ \ \\`\\ \ \//\ \L\.\_  
  \ \_\ \_\ \_\ \_\\ \__/.\_\     
   \/_/\/_/\/_/\/_/ \/__/\/_/
                             
EOF

# Usage
function usage() {
  cat << EOS >&2        
ikra ${VERSION} -RNAseq pipeline centered on Salmon-
Usage: ${PROGNAME} experiment_table.csv [--test --fastq, --help, --without-docker, --udocker, --protein-coding] [--threads [VALUE]][--output [VALUE]][--suffix_PE_1 [VALUE]][--suffix_PE_2 [VALUE]][--transcript [VALUE]][--genome [VALUE]]
  args
    1.experiment matrix(csv)
    2.reference(human or mouse)
Options:
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

Citation :
Hiraoka, Y., Yamada, K., Yamasaki, R., Kawasaki, Y., Kitabatake, R., Matsumoto, Y., Ishikawa, K., Umezu, Y., Hirose, H., & Yasumizu, Y. (2021). ikra v2.0: RNAseq pipeline centered on Salmon. https://doi.org/10.5281/zenodo.4718200

Github repo : https://github.com/yyoshiaki/ikra
EOS
  exit 1
}