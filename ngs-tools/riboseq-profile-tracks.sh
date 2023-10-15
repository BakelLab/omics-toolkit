#!/bin/sh

# 04.12.2019 09:59:53 EST
# Harm van Bakel <hvbakel@gmail.com>

SAMPLE=""
while getopts "i:" opt; do
   case $opt in
   
   i)
      SAMPLE="$OPTARG"
      ;;
   
   *)
      echo "Incorrect options provided"
      exit 1
      ;;

   esac
done

# Check arguments and produce a help message
if [ -z "$SAMPLE" ];
then
  cat << EOF

   Usage: analysis_p-site-distribution.sh -i <sample-name>

   Arguments:
    -i <string>
      Name of sample
    -help
      This help message

EOF
  exit 0
fi

# Load required modules
module load R

# Set global aliases
alias awkt="awk -F '\t' -v OFS='\t'"
alias sortt="sort -t $'\t'"

# Check if required files exist
if [ ! -e "${SAMPLE}.fastq.gz" ]; then
   echo "Required input file '${SAMPLE}.fastq.gz' does not exist"
   exit 0;
fi 
if [ ! -e "${SAMPLE}_genome.init.pos.bw" ]; then
   echo "Required input file '${SAMPLE}.fastq.gz' does not exist"
   exit 0;
fi 

#######################################
# EXTRACT CS EVENTS FROM RIBOSEQ DATA #
#######################################

# Grab CS events for each segment (with 7 nt of segment-specific sequence) and remove any canonical 5' cRNA hit
gunzip -c ${SAMPLE}.fastq.gz | grep -P 'GC[AG]AAAGCAGGTCA' | perl -pe 's/AAAAAAAA.*$//' | grep -vP '^AGC[AG]AAA' | awkt '{print $1, "seg1"}' >  ${SAMPLE}_CS-events.txt
gunzip -c ${SAMPLE}.fastq.gz | grep -P 'GC[AG]AAAGCAGGCAA' | perl -pe 's/AAAAAAAA.*$//' | grep -vP '^AGC[AG]AAA' | awkt '{print $1, "seg2"}' >> ${SAMPLE}_CS-events.txt
gunzip -c ${SAMPLE}.fastq.gz | grep -P 'GC[AG]AAAGCAGGTAC' | perl -pe 's/AAAAAAAA.*$//' | grep -vP '^AGC[AG]AAA' | awkt '{print $1, "seg3"}' >> ${SAMPLE}_CS-events.txt
gunzip -c ${SAMPLE}.fastq.gz | grep -P 'GC[AG]AAAGCAGGGGA' | perl -pe 's/AAAAAAAA.*$//' | grep -vP '^AGC[AG]AAA' | awkt '{print $1, "seg4"}' >> ${SAMPLE}_CS-events.txt
gunzip -c ${SAMPLE}.fastq.gz | grep -P 'GC[AG]AAAGCAGGGTA' | perl -pe 's/AAAAAAAA.*$//' | grep -vP '^AGC[AG]AAA' | awkt '{print $1, "seg5"}' >> ${SAMPLE}_CS-events.txt
gunzip -c ${SAMPLE}.fastq.gz | grep -P 'GC[AG]AAAGCAGGGGT' | perl -pe 's/AAAAAAAA.*$//' | grep -vP '^AGC[AG]AAA' | awkt '{print $1, "seg6"}' >> ${SAMPLE}_CS-events.txt
gunzip -c ${SAMPLE}.fastq.gz | grep -P 'GC[AG]AAAGCAGGTAG' | perl -pe 's/AAAAAAAA.*$//' | grep -vP '^AGC[AG]AAA' | awkt '{print $1, "seg7"}' >> ${SAMPLE}_CS-events.txt
gunzip -c ${SAMPLE}.fastq.gz | grep -P 'GC[AG]AAAGCAGGGTG' | perl -pe 's/AAAAAAAA.*$//' | grep -vP '^AGC[AG]AAA' | awkt '{print $1, "seg8"}' >> ${SAMPLE}_CS-events.txt

# Get the lengths of the host portions for all CS reads
cat ${SAMPLE}_CS-events.txt | perl -pe 's/.GC[AG]AAAGCA.*//' | awk '{print "-" length($1)}' > ${SAMPLE}_CS-events.tmp1
cut -f 2 ${SAMPLE}_CS-events.txt > ${SAMPLE}_CS-events.tmp2
paste ${SAMPLE}_CS-events.tmp1 ${SAMPLE}_CS-events.tmp2 > ${SAMPLE}_CS-events_lengths.txt
rm -f ${SAMPLE}_CS-events.tmp1 ${SAMPLE}_CS-events.tmp2

# Get the lengths of the host portions for ATG events only
cat ${SAMPLE}_CS-events.txt | perl -pe 's/GC[AG]AAAGCA.*//' | grep ATG | awk '{print "-" length($1)}' > ${SAMPLE}_CS-events_AUG.tmp1
cat ${SAMPLE}_CS-events.txt | perl -pe 's/GC[AG]AAAGCA.*\t/\t/' | grep ATG | cut -f 2 > ${SAMPLE}_CS-events_AUG.tmp2
paste ${SAMPLE}_CS-events_AUG.tmp1 ${SAMPLE}_CS-events_AUG.tmp2 > ${SAMPLE}_CS-events_lengths_AUG.txt
rm -f ${SAMPLE}_CS-events_AUG.tmp1 ${SAMPLE}_CS-events_AUG.tmp2

# Get histogram of lengths
awkt '{out[$2 "_" $1]++} END{for(var in out){print var, out[var]}}' ${SAMPLE}_CS-events_lengths.txt     | perl -pe 's/_/\t/' | sortt -k1,1 -k2,2n > ${SAMPLE}_CS-events_lengths_hist.txt 
awkt '{out[$2 "_" $1]++} END{for(var in out){print var, out[var]}}' ${SAMPLE}_CS-events_lengths_AUG.txt | perl -pe 's/_/\t/' | sortt -k1,1 -k2,2n > ${SAMPLE}_CS-events_lengths_AUG_hist.txt

##############################################
# GET OTHER VIRAL SEGMENT MAPPED READ STARTS #
##############################################

# Convert mapped positive-sense reads to histogram; restrict to just the first 80 nt of the viral segments
bigWigToBedGraph ${SAMPLE}_genome.init.pos.bw ${SAMPLE}_SG-events.tmp
grep seg ${SAMPLE}_SG-events.tmp | awkt '$2<=80 {print $1, $2,$4}' > ${SAMPLE}_SG-events_hist.txt
rm -f ${SAMPLE}_SG-events.tmp

####################
# MAKE FACET PLOTS #
####################

# Plot event lengths in a faceted ggplot
cat  ${SAMPLE}_SG-events_hist.txt ${SAMPLE}_CS-events_lengths_hist.txt | sortt -k1,1 -k2,2n > ${SAMPLE}_riboseq-histogram.txt
add-header ${SAMPLE}_riboseq-histogram.txt Segment Position Count
plot-riboseq-histograms.R -i ${SAMPLE}_riboseq-histogram.txt -o ${SAMPLE}

# Get per-segment CS event counts with/without AUGs
awkt '{out[$2]++} END{for(var in out){print var, out[var]}}' ${SAMPLE}_CS-events_lengths.txt | sortt -k1,1     > ${SAMPLE}_count_IAV-CS
awkt '{out[$2]++} END{for(var in out){print var, out[var]}}' ${SAMPLE}_CS-events_lengths_AUG.txt | sortt -k1,1 > ${SAMPLE}_count_IAV-CS-AUG
awkt '{out[$1] = out[$1] + $3} END{for(var in out){print var, out[var]}}' ${SAMPLE}_riboseq-histogram.txt | grep -v Segment | sortt -k1,1 > ${SAMPLE}_count_IAV-total
multi-join-by-ids -e 0 ${SAMPLE}_count_IAV-CS-AUG ${SAMPLE}_count_IAV-CS ${SAMPLE}_count_IAV-total | sortt -k1,1 > ${SAMPLE}_riboseq-summary.txt


#########################################
# COUNT INFLUENZA AND HOST MAPPED READS #
#########################################

bigWigToBedGraph ${SAMPLE}_genome.init.pos.bw ${SAMPLE}_SG-events.tmp
iavscount=`awkt '$1~/^seg/ {out+=$4} END{print out}' ${SAMPLE}_SG-events.tmp`
iavccount=`awkt '$1~/^seg/ {out+=$3} END{print out}' ${SAMPLE}_CS-events_lengths_hist.txt`
iavtcount=$(($iavscount + $iavccount))
hostcount=`awkt '$1!~/^seg/ {out+=$4} END{print out}' ${SAMPLE}_SG-events.tmp`
sumcount=$(($iavtcount + $hostcount))
echo -e "$SAMPLE\t$iavtcount\t$hostcount\t$sumcount" > ${SAMPLE}_host-virus.total.counts


##########################################
# GET MAPPED READ LENGTHS FOR VIRUS/HOST #
##########################################

# Get histograms
bamToBed -bed12 -i ${SAMPLE}_transcriptome.bam | awkt '$10==1 && $1~"^ENST" && $4!=last && $6=="+" {out[$11]++; last=$4} END{for(var in out){print var, out[var]}}' | sortt -k1,1 -k2,2n > ${SAMPLE}_mapped-lengths_host.txt
bamToBed -bed12 -i ${SAMPLE}_transcriptome.bam | awkt '$10==1 && $1~"_pos$" && $4!=last && $6=="+" {out[$11]++; last=$4} END{for(var in out){print var, out[var]}}' | sortt -k1,1 -k2,2n > ${SAMPLE}_mapped-lengths_virus.txt

############################
# GET P-SITE OFFSET TRACKS #
############################

