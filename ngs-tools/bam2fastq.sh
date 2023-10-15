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

# Figure out how many fastq files were aligned (i.e. distinguish between paired and single-end
fastq_count=`samtools view -H ${bamfile} | grep '^@PG'| head -n1 | grep -o '\.fastq' | wc -l`

# Set parameters for terashuf
export TMPDIR=/sc/arion/scratch/$USER/tmp/
export MEMORY=20
echo "TMPDIR: $TMPDIR" >  ${outpath}/${outprefix}_bam2fastq.log
echo "MEMORY: $MEMORY" >> ${outpath}/${outprefix}_bam2fastq.log

# Do conversion for single or paired-end format
if [[ $fastq_count = 2 ]]
then
   java -jar ${PICARD} SamToFastq \
      INPUT=${bamfile} \
      FASTQ=${outpath}/${outprefix}_1.fastq.gz \
      SECOND_END_FASTQ=${outpath}/${outprefix}_2.fastq.gz \
      UNPAIRED_FASTQ=${outpath}/${outprefix}_up.fastq.gz \
      VALIDATION_STRINGENCY=SILENT \
      >> ${outpath}/${outprefix}_bam2fastq.log 2>&1

   paste <(zcat ${outpath}/${outprefix}_1.fastq.gz) <(zcat ${outpath}/${outprefix}_2.fastq.gz) | paste - - - - | terashuf | awk -v outa=${outpath}/${outprefix}_1.fastq.shuf -v outb=${outpath}/${outprefix}_2.fastq.shuf -F'\t' '{OFS="\n"; print $1,$3,$5,$7 | "gzip > "outa".gz"; print $2,$4,$6,$8 | "gzip > "outb".gz"}'
   rm -f ${outpath}/${outprefix}_1.fastq.gz ${outpath}/${outprefix}_2.fastq.gz
else
   java -jar ${PICARD} SamToFastq \
      INPUT=${bamfile} \
      FASTQ=${outpath}/${outprefix}.fastq.gz \
      VALIDATION_STRINGENCY=SILENT \
      >> ${outpath}/${outprefix}_bam2fastq.log 2>&1

   gunzip -c ${outpath}/${outprefix}.fastq.gz | terashuf | gzip > ${outpath}/${outprefix}.fastq.shuf.gz
   rm -f ${outprefix}.fastq.gz
fi
