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

# Check for picard
if [[ -z "$PICARD" ]]
then
   if ! command -v picard &> /dev/null
   then
      if command -v module &> /dev/null
      then
         module load picard && echo "Loaded picard module" || echo "Error: failed to load picard module"
      else
         echo "Error: picard is not installed"
         exit 1
      fi
   else
      PICARD="picard"
   fi
fi

# Do conversion for single or paired-end format
java -jar ${PICARD} SamToFastq \
   INPUT=${bamfile} \
   FASTQ=${outpath}/${outprefix}.fastq.gz \
   VALIDATION_STRINGENCY=SILENT \
   >> ${outpath}/${outprefix}_bam2fastq.log 2>&1
