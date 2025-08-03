#!/bin/bash

# --- Configuration ---
# This script converts a BAM file to shuffled FASTQ, auto-detecting
# if the reads are single-end or paired-end. It uses terashuf for
# out-of-core shuffling and pigz for multi-threaded compression.

set -e # Exit immediately if a command exits with a non-zero status.

# --- USAGE ---
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 /path/to/your.bam output_prefix"
    echo "Example: $0 my_data.bam my_shuffled_reads"
    exit 1
fi

# --- Variables ---
BAM_FILE="$1"
OUT_PREFIX="$2"
THREADS=12
RANDOM_SEED=100 # Use a fixed seed for reproducibility

# --- Dependency Checks ---
if ! command -v terashuf &> /dev/null; then
    echo "Error: terashuf is not installed" >&2
    exit 1
fi
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools is not installed" >&2
    exit 1
fi
if ! command -v pigz &> /dev/null; then
    echo "Error: pigz is not installed" >&2
    exit 1
fi

# --- Paired-End Detection ---
echo "🔎 Checking for paired-end reads..."
PAIRED_COUNT=$(samtools view -c -f 1 -F 4 "$BAM_FILE")

# --- Setup for Terashuf ---
export TMPDIR=${SCRATCH:-/tmp}/$USER/tmp/
export MEMORY=20
mkdir -p "$TMPDIR"
echo "TMPDIR: $TMPDIR" >  "${OUT_PREFIX}_bam2fastq.log"
echo "MEMORY: $MEMORY" >> "${OUT_PREFIX}_bam2fastq.log"

# --- Main Logic ---
if [[ $PAIRED_COUNT -gt 0 ]]; then
    ### PAIRED-END WORKFLOW ###
    echo "✅ Paired-end reads detected. Starting fully streaming conversion..."
    
    samtools collate -@ $THREADS -O "$BAM_FILE" | \
      samtools fastq -@ $THREADS -n - | \
      paste - - - - - - - - | \
      terashuf | \
      awk -v outa="${OUT_PREFIX}_1.fastq.shuf" -v outb="${OUT_PREFIX}_2.fastq.shuf" -v p="$THREADS" -F'\t' '{OFS="\n"; print $1,$2,$3,$4 | "pigz -p "p" > "outa".gz"; print $5,$6,$7,$8 | "pigz -p "p" > "outb".gz"}'

else
    ### SINGLE-END WORKFLOW ###
    echo "✅ No paired reads found. Running single-end workflow..."
    
    samtools fastq -@ $THREADS "$BAM_FILE" | \
      terashuf | \
      pigz -p "$THREADS" > "${OUT_PREFIX}.fastq.shuf.gz"
fi

echo "🎉 Done."
