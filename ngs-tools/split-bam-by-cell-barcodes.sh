#!/bin/sh

# 28.02.2021 10:59:45 EST
# Harm van Bakel <hvbakel@gmail.com>

# Globals
bamfile=$1
barcodefile=$2
outprefix=$3

# Help message
if [[ -z "$1" ]]
then
   echo -e "\nUsage: split-bam-by-cell-barcodes.sh <bam> <barcodefile> <outprefix>\n"
   exit 0
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

# Set TMPDIR
export TMPDIR=/sc/arion/scratch/$USER/tmp/


