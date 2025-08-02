#!/bin/bash

# 01.09.2018 14:10:58 EDT
# Harm van Bakel <hvbakel@gmail.com>

# Globals
bamfile=$1
outpath=${2:-./}
outprefix=`basename $bamfile .bam`

if [[ -z "$1" ]]
then
   echo -e "\nUsage: bam2fastq.sh <bam> [output-path]\n"
   exit 0
fi

# Check for terashuf
if ! command -v terashuf &> /dev/null
then
   echo "Error: terashuf is not installed"
   exit 1
fi

# Check for samtools
if ! command -v samtools &> /dev/null
then
   if command -v module &> /dev/null
   then
      module load samtools && echo "Loaded samtools module" || echo "Error: failed to load samtools module"
   else
      echo "Error: samtools is not installed"
      exit 1
   fi
fi

# --- Paired-End Detection ---
# A more robust check for paired-end data is to count reads with the 'paired' flag.
# We also skip unmapped reads, as they might not have the flag set correctly.
echo "🔎 Checking for paired-end reads..."
fastq_count=$(samtools view -c -f 1 -F 4 "$bamfile")

# --- Setup for Terashuf ---
export TMPDIR=${SCRATCH:-/tmp}/$USER/tmp/
export MEMORY=20
mkdir -p "$TMPDIR"
echo "TMPDIR: $TMPDIR" >  "${outpath}/${outprefix}_bam2fastq.log"
echo "MEMORY: $MEMORY" >> "${outpath}/${outprefix}_bam2fastq.log"

# Do conversion for single or paired-end format
if [[ $fastq_count -gt 0 ]]; then
    echo "✅ Paired-end reads detected. Starting fully streaming conversion..."
    
    # This single pipeline does everything in memory and on-the-fly.
    samtools fastq -@ 8 -n "$bamfile" | \
    paste - - - - - - - - | \
    terashuf | \
    awk -v outa="${outpath}/${outprefix}_1.fastq.shuf" -v outb="${outpath}/${outprefix}_2.fastq.shuf" -F'\t' '{OFS="\n"; print $1,$2,$3,$4 | "gzip > "outa".gz"; print $5,$6,$7,$8 | "gzip > "outb".gz"}'

else
    # The single-end logic is already fully streaming and efficient.
    echo "✅ Single-end reads detected. Starting conversion..."
    samtools fastq -@ 8 "$bamfile" | terashuf | gzip > "${outpath}/${outprefix}.fastq.shuf.gz"
fi
