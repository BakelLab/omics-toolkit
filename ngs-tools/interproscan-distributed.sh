#!/bin/sh

# 08.07.2015 13:24:02 EDT
# Harm van Bakel <hvbakel@gmail.com>

# Check arguments
if [ $# -lt 3 ]
then
  echo "Usage: `basename $0` <aa-fasta-file> <seqs-per-batch> <minerva-allocation> <walltime> <extra-interpro-args>"
  exit 0
fi

# Set path to parent dir for git repositories
GIT_REPODIR="${GIT_REPODIR:=/hpc/users/vanbah01/opt}"
export GIT_REPODIR

# Set paths to minerva-queue-lsf repository
if [ -d "$GIT_REPODIR/minerva-queue-lsf/" ]; then
   export PATH="$GIT_REPODIR/minerva-queue-lsf/bin:$PATH"
else
   echo "Could not find the minerva-queue-lsf repository in '$GIT_REPODIR'. Check GIT_REPODIR environment variable."
   exit 1
fi

# Set paths to ngs-tools repository
if [ -d "$GIT_REPODIR/ngs-tools/" ]; then
   export PATH="$GIT_REPODIR/ngs-tools/bin:$PATH"
   export PATH="$GIT_REPODIR/ngs-tools/pipelines/:$PATH"
else
   echo "Could not find the ngs-tools repository in '$GIT_REPODIR'. Check GIT_REPODIR environment variable."
   exit 1
fi

# Set paths to igb-tools repository
if [ -d "$GIT_REPODIR/igb-tools/" ]; then
   export PATH="$GIT_REPODIR/igb-tools/bin:$PATH"
else
   echo "Could not find the igb-tools repository in $GIT_REPODIR. Check GIT_REPODIR environment variable."
   exit 1
fi

# Set path to utility repository
if [ -d "$GIT_REPODIR/utility/" ]; then
   export PATH="$GIT_REPODIR/utility/bin:$PATH"
   export PERL5LIB="$GIT_REPODIR/utility/perl5"
else
   echo "Could not find the igb-tools repository in $GIT_REPODIR. Check GIT_REPODIR environment variable."
   exit 1
fi

# Reset module system for node architecture
module purge all

# Check if we have a module system
if type "module" >/dev/null 2>&1;
then
   # Load modules
   module load java/11.0.2
   module load python
fi

# Set output filename
OUT=`basename $1 .fa`
OUT=`basename $OUT .faa`
OUT=`basename $OUT .fasta`

# Generate a random name for the jobs
JOBNAME=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 5 | head -n 1)

# Split input fasta according to the number of seqs per file
seqcount=`grep -c '>' $1`
let batchcount=$seqcount/$2
echo "Submitting $batchcount jobs"
$GIT_REPODIR/ngs-tools/bin/fasta-splitter-random.sh $1 $batchcount ${JOBNAME}

# Submit the jobs to run
COUNTER=0
JOBCONDITION=""
JOBERRORCHECK=""
export TMPDIR="/sc/arion/scratch/${USER}/"
export TMP_DIR="/sc/arion/scratch/${USER}/"

for i in ${JOBNAME}_*.fa
do
   name=`basename $i .fa`
   $GIT_REPODIR/minerva-queue-lsf/bin/submitjob $4 -c 12 -A $3 -q premium -J ${JOBNAME}${COUNTER}_ips \
     /sc/arion/projects/bakellab/reference-databases/interproscan/interproscan-5.68-100.0/interproscan.sh -T /sc/arion/scratch/${USER}/ -i $i -goterms ${@:5} \> $name.out

   # Increment job condition and counter
   if [ -z "$JOBCONDITION" ];
   then
      JOBCONDITION="ended(\"${JOBNAME}${COUNTER}_ips\")"
      JOBERRORCHECK="exit(\"${JOBNAME}${COUNTER}_ips\",!=0)"
   else
      JOBCONDITION="$JOBCONDITION && ended(\"${JOBNAME}${COUNTER}_ips\")"
      JOBERRORCHECK="$JOBERRORCHECK && exit(\"${JOBNAME}${COUNTER}_ips\",!=0)"
   fi
   COUNTER=$((COUNTER + 1))
done

# Submit the final merge job with a dependency on the jobs that were just submitted
$GIT_REPODIR/minerva-queue-lsf/bin/submitjob 1 -c 1 -A $3 -q premium -J gather_${JOBNAME} -w "${JOBCONDITION}" \
   cat ${JOBNAME}_\*.fa.tsv \> ${OUT}.tsv \; \
   cat ${JOBNAME}_\*.fa.gff3 \> ${OUT}.gff3 \; \
   cat ${JOBNAME}_\*.out \> ${OUT}.log \; \
   iprscan2topgoinput.pl ${OUT}.tsv \> ${OUT}.topGO \; \
   rm -rf ${JOBNAME}_\*

# Make sure we know if any of the individual jobs failed with an error by submitting another job that only runs if there was an error in the parent
$GIT_REPODIR/minerva-queue-lsf/bin/submitjob 1 -c 1 -A $3 -q premium -J check_${JOBNAME} -w "${JOBERRORCHECK}" \
   echo "One or more jobs failed with an error, please rerun the interpro analysis" \> ${OUT}.ERROR
