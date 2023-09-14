#!/bin/bash
# split a fasta file into multiple files, one file per sequence in the input fasta file
# use the first field of the sequence id (white-space separated) as name for the output file
# input argument 1 is gzipped fasta file, 
# input argument 2 is output directory
infile=$(realpath $1)
outdir=$2

# Version 1: one sequence per file
mkdir -p ${outdir} \
  && cd $(realpath ${outdir}) \
  && gunzip -c ${infile} \
  | awk '/^>/{split(substr($0,2),a," "); filename=a[1] ".fa.gz"}; {print | "gzip > " filename; close(filename)}' \
  && cd -

exit 0

# Version 2: one chromosome "NC" per file, the rest into one file DOES NOT WORK YET !
mkdir -p ${outdir} \
  && cd $(realpath ${outdir}) \
  && rm -f nonchromosomal.fa.gz \
  && gunzip -c ${infile} \
  | awk '/^>NC_/{split(substr($0,2),a," "); filename=a[1] ".fa.gz"}; {print | "gzip > " filename; close(filename)}; \
         /^>NW_/{split(substr($0,2),a," "); filename="nonchromosomal.fa.gz"}; {print | "gzip >> " filename; close(filename)}' \
  && cd -
  
exit 0

