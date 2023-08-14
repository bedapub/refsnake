#!/bin/bash

# ----------------------------
# Make rRNA inverals file
# ----------------------------

if [[ $# -ne 3 ]]; then
  echo "3 input arguments required, in this order:"
  echo "  1. fasta/genome.chrom.sizes"
  echo "  2. gff3/refseq/refseq.gff3"
  echo "  3. gff3/ensembl/ensembl.gff3"
  echo ""
  echo "Results go to stdout"
  exit 0
fi

input_sizes=$(realpath $1)
input_refseq=$(realpath $2)
input_ensembl=$(realpath $3)

if [ ! -e ${input_sizes} ]; then echo "Input file ${input_size} does not exist. Abort!" && exit 1; fi
if [ ! -e ${input_refseq} ]; then echo "Input file ${input_refseq} does not exist. Abort!" && exit 1; fi
if [ ! -e ${input_ensembl} ]; then echo "Input file ${input_ensembl} does not exist. Abort!" && exit 1; fi


file_size=$(stat -c "%s" ${input_sizes})
if [ ${file_size} -eq 0 ]; then echo "Input file ${input_size} is empty. Abort!" && exit 1; fi

tmpfile=$(mktemp)
output_tmp=${tmpfile}.tmp
output_header=${tmpfile}.header
output_refseq=${tmpfile}.refseq
output_ensembl=${tmpfile}.ensembl

# Output chromosome size as SAM header
perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]"' ${input_sizes} > ${output_header}
  
# Output intervals for rRNA transcripts
# From Ensembl
file_size=$(stat -c "%s" ${input_ensembl})
if [ ${file_size} -eq 0 ]; then
  touch ${output_ensembl}
else
  chroms="$(grep -v '^#' ${input_ensembl} | cut -f1 | sort | uniq)"
  rm -f ${output_ensembl}
  for i in ${chroms}; do  
    grep -w ${i} ${input_ensembl} | awk 'BEGIN{FS="\t"}{if ($3=="rRNA"){print $0}}' | cut -f1,4,5,7,9 \
      | perl -lane '/Name=([^"]+)/ or die "no transcript_id on $."; print join "\t", (@F[0,1,2,3], $1)' \
      | cut -f1 -d\; | sort -k1V -k2n -k3n >> ${output_ensembl}
  done
fi
  
# From RefSeq
file_size=$(stat -c "%s" ${input_refseq})
if [ ${file_size} -eq 0 ]; then
  touch ${output_refseq}
else
  chrom="$(grep -v '^#' ${input_refseq} | cut -f1 | sort | uniq)"
  rm -f ${output_refseq}
  for i in ${chroms}; do
    grep -w ${i} ${input_refseq} | awk 'BEGIN{FS="\t"}{if ($3=="rRNA"){print $0}}' | cut -f1,4,5,7,9 \
      | perl -lane '/gene=([^"]+)/ or die "no transcript_id on $."; print join "\t", (@F[0,1,2,3], $1)' \
      | cut -f1 -d\; | sort -k1V -k2n -k3n >> ${output_refseq}
  done
fi
  
cat ${output_ensembl} ${output_refseq} | sort -k1V -k2n -k3n | uniq > ${output_tmp}
cat ${output_header} > ${tmpfile}
awk 'BEGIN{FS="\t"}{s=sprintf("%s\t%s\t%s\t%s",$1,$2,$3,$4); arr[s]=sprintf("%s,%s",arr[s],$5)} END \
  {for (i in arr) {printf("%s\t%s\n",i,arr[i])}}' ${output_tmp} | sed 's/\t,/\t/g' \
  | sort -k1V -k2n -k3n >> ${tmpfile}

rm -f ${output_tmp} ${output_refseq} ${output_ensembl} ${output_header}
cat ${tmpfile}


exit 0

# Perl code translated into Python by ChatGPT: 
# Need to be checked first:
import sys
import re

for line in sys.stdin:
    line = line.strip()
    fields = line.split('\t')
    match = re.search(r'Name=([^"]+)', line)
    if not match:
        sys.stderr.write(f'no transcript_id on line {line}\n')
        continue
    name_value = match.group(1)
    output = '\t'.join(fields[:4] + [name_value])
    print(output)