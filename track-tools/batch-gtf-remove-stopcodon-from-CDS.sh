#!/bin/bash

# 14.06.2020 11:55:14 EDT
# Harm van Bakel <hvbakel@gmail.com>

# Check usage
if (( $# != 2 )); then
    echo -e "\n  Usage: batch-gtf-remove-stopcodon-from-CDS.sh <isoseq.gtf> <output.gtf>\n"
    exit;
fi

# Load anaconda module
# conda activate IsoformSwitchAnalyzeR

# Generate a random name for the jobs
JOBNAME=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 5 | head -n 1)

# Split input file by chromosome
chr_list=`cut -f 1 $1 | uniq | sort | uniq | grep chr`
for i in ${chr_list}
do
   awk -F '\t' -v OFS='\t' -v chr=$i '$1==chr' $1 > input_${JOBNAME}_${i}.gtf
done

# Run mikado for 4 genomes at a time in parallel
parallel --jobs 24 "gtf-remove-stopcodon-from-CDS.R -i input_${JOBNAME}_{}.gtf -o output_${JOBNAME}_{}.gtf 2> output_${JOBNAME}_{}.log" ::: ${chr_list}

# Combine final outputs
cat output_${JOBNAME}_*.gtf > $2
