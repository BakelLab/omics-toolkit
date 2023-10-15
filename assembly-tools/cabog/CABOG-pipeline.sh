#!/bin/sh

# 01.12.2010 12:59:35 EST
# Harm van Bakel <hvbakel@gmail.com>


#---------------------------------------
# 0-gkpStore & 0-mercounts
#---------------------------------------

# Submit cabog run (initial store building)
submitjob 24 -c 16 -m 120 runCA -d $PWD/assembly -p test -s runca.spec stopAfter=initialStoreBuilding

# Now run initialTrim and meryl (note that it is possible to submit this job with a dependency on the previous one)
submitjob 48 -c 16 -m 120 runCA -d $PWD/assembly -p test -s runca.spec stopAfter=overlapper

# It is very likely that cabog sets the kmer maxCount threshold way too high. Check the
# meryl dir and look at the value in the .estMerThresh.out file. Run the following command 
# (with modified paths etc.) to get a fasta file of kmers to skip at other thresholds.
# Run this in the meryl folder
meryl -Dt -n 250 -s test-C-ms22-cm0 > test.nmers.obt.fasta.n600 2> test.nmers.obt.fasta.n600.err


#---------------------------------------
# 0-overlaptrim-overlap & 0-overlaptrim
#---------------------------------------

# After meryl and the initial trimming is done, the overlap trimming stage begins
# Since cabog doesn't integrate with the moab system, we'll need to manually submit the overlapper
# jobs. First make sure to start runCA to generate the required files in the 0-overlaptrim-overlap
# dir: 'overlap.sh', 'ovlbat', 'ovljob', and 'ovlopt'. Run the following command to generate the 
# required files:
submitjob 48 -c 16 -m 120 runCA -d $PWD/assembly -p test -s runca.spec useGrid=1 scriptOnGrid=0 stopAfter=overlapper

# Once the required files are there, change into the '0-overlaptrim-overlap' and run the following
# command to run the first overlapper stage. This script submits the overlapper jobs in such a way
# that only 80 (default, this can be changed) overlapper jobs run at any given time, so that the
# IO doesn't swamp the cluster disks. NOTE: you may need to edit the overlap.sh script to make sure
# that it points to the correct kmer exclusion file in the 0-mercounts folder.
~/bin/submit-overlapper-jobs.pl

# Once all overlapper jobs finish, it's time to build the first overlapstore and run several associated
# data processing steps (overlap-based trimming, and detecting chimeric and spur reads). 
# NOTE: submit the second job with dependency on completion of the first job, that way they will run 
# sequentially unattended.
submitjob 48 -c 16 -m 120 runCA -d $PWD/assembly -p test -s runca.spec stopBefore=chimeraDetection
submitjob 48 -c 16 -m 120 -W depend=afterany:LASTJOB runCA -d $PWD/assembly -p test -s runca.spec stopAfter=OBT

#--------------------------
# 1-overlapper
#--------------------------
# At this point the first overlap store will be ready and the gatekeeper store has been updated with new
# clear ranges (i.e. the thresholds on what constitutes reliable sequence). Now we need to rerun overlapper
# in exactly the same way as in the first stage, with the difference that the updated clear ranges are used
# this time. Since some reads will be excluded in the previous stage, the overlap store in this stage will
# be significantly smaller than in the first run. First generate the required files for job submission in
# the 1-overlapper dir ('overlap.sh', 'ovlbat', 'ovljob', and 'ovlopt')
submitjob 48 -c 16 -m 120 runCA -d $PWD/assembly -p test -s runca.spec useGrid=1 scriptOnGrid=0

# Once the required files are there, change into the '1-overlapper' dir and run the following
# command to run the second overlapper stage. This script submits the overlapper jobs in such a way
# that only 80 (default, this can be changed) overlapper jobs run at any given time, so that the
# IO doesn't swamp the cluster disks. NOTE: you may need to once again edit the overlap.sh script to 
# make sure that it points to the correct kmer exclusion file in the 0-mercounts folder.
~/bin/submit-overlapper-jobs.pl

# Wait for all the overlapper jobs to finish and then proceed to building the second overlapstore:
submitjob 48 -c 16 -m 120 runCA -d $PWD/assembly -p test -s runca.spec stopAfter=overlapper

#--------------------------
# 3-overlapcorrection
#--------------------------

# Examine fragment multialigns for errors. This stage can once again be parallellized, so we first generate
# the run script. Note that the number of parallel jobs should be kept small because of the amount of IO
runCA -d $PWD/assembly -p test -s runca.spec useGrid=1 scriptOnGrid=0

# The previous command will indicating how many jobs are needed. The fragment correction is optimized to run
# on high memory nodes, since each batch that is processed requires reading the entire gatekeeper store.
# A batch of 50 million illumina reads requires approximately 110 Gb of ram. Since the number of jobs is limited,
# it's possible to simply submit by hand. This stage takes 18-36 hours, depending on largemem node availability
cd $PWD/assembly/3-overlapcorrection
for i in `seq 1 8`
do
   submitjob 6 -c 16 -m 120 ./frgcorr.sh $i
done

# Apply error corrections to overlaps. This stage starts with the concatenation of all the frgcorr files from
# the previous stage, so we'll need to submit this to a node. After that, there are a bunch of parallel jobs again.
submitjob 12 -c 16 -m 120 runCA -d $PWD/assembly -p test -s runca.spec stopBefore=unitigger

# Now we submit the jobs. These jobs are much less IO intensive, so we can run on a bunch of standard nodes.
# For 2.5M reads, each job takes about 45 mins. At 80 concurrent jobs, the whole thing should finish in 4 hours
~/bin/submit-ovlcorr-jobs.pl

# This stage ends with a single-threaded process to update the overlap store with corrected error rates.
# It takes about 38 hours to update a 1Tb overlapstore.
submitjob 48 -c 16 -m 120 runCA -d $PWD/assembly -p test -s runca.spec stopBefore=unitigger

###################################################
# ===>>>  MOVE SEQUENCE STORES TO BEHEMOTH <<<=== #
###################################################

#--------------------------
# 4-unitigger
#--------------------------

# Unitigger and all subsequent stages are run on the 'behemoth' largemem server. Gzip the sequence stores on
# scinet and then copy them over to behemoth using rsync on datamover1.
# The unitig stage builds the tigstore. If unitigger gets killed prematurely, it's possible to pick up 
# where it left off
runCA -d $PWD/assembly -p test -s runca.spec > runca.1.err 2>&1


#--------------------------
# 5-consensus
#--------------------------

# This stage computes the unitig consensus sequences. Control the number of jobs on behemoth using the 
# 'cnsConcurrency' parameter. This stage should start directly after unitigging, but if there was a problem
# during unitigging, simply restart with:
runCA -d $PWD/assembly -p test -s runca.spec > runca.2.err 2>&1

# There will always be some segments that fail during consensus sequence building
# These can be fixed with the following script
~/hugheslab/svncheckout/assembly/trunk/cabog/fix-unitigs.pl

#--------------------------
# 7-CGW
#--------------------------

# Once all unitig consensus sequences have been collected, it's time for the scaffolding phase
# Again, it's simply a matter of starting runCA after all the unitigs have been fixed
runCA -d $PWD/assembly -p test -s runca.spec > runca.3.err 2>&1


#--------------------------
# 8-consensus
#--------------------------

# Basically the same as stage 5, but now for scaffolds rather than contigs.
# Note that the memory usage of each consensus sequence process is much larger than during the
# unitigging stage, so make sure not to start too many jobs in parallel (cnsConcurrency)
# Consensus should run automatically after the scaffolding stage, but it's possible to restart if necessary
runCA -d $PWD/assembly -p test -s runca.spec > runca.4.err 2>&1

#--------------------------
# 9-terminator
#--------------------------

# The final stage collects all the unitig and scaffolding info and writes it to disk together with QC info
# Again, this should run automatically after 8-consensus, however, it's possible to restart if there
# was a problem during scaffolding.
runCA -d $PWD/assembly -p test -s runca.spec > runca.5.err 2>&1
