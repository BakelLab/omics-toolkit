#!/bin/sh

# 20.02.2017 17:31:11 EST
# Harm van Bakel <hvbakel@gmail.com>


grep -vP '^@'  ~/projects/reference-databases/RNASeq-pipeline/hg19-v19-ERCC/SeqCap_lncRNA_hg19_capture_targets.list > SeqCap_lncRNA_hg19_capture_targets.bed
alias awkt="awk -F '\t' -v OFS='\t'"


for i in *.sam
do
   name=`basename $i .sam`
   samtools view -b ${i} > ${name}.tmp
   samtools sort ${name}.tmp ${name}
   samtools index ${name}.bam
   rm -f ${name}.tmp
   intersectBed -u -a ${name}.bam -b SeqCap_lncRNA_hg19_capture_targets.bed > ${name}_in.bam
   samtools view ${name}_in.bam | cut -f 1 | perl -pe 's/.*\/f//; s/\/.*//; s/p/\t/' | awkt '{out+=$1+$2} END{print out}' > ${name}_in.count
   samtools view ${name}.bam | cut -f 1 | perl -pe 's/.*\/f//; s/\/.*//; s/p/\t/' | awkt '{out+=$1+$2} END{print out}' > ${name}.count
done
