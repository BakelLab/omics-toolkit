#!/usr/bin/env perl

# 09.09.2018 11:43:30 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(cwd);
use IPC::Cmd qw[can_run run];

# GET PARAMETERS
my $sHelp          = 0;
my $sInput         = "";
my $sGenome        = "";
my $sStrand        = "opposite";
my $sTask          = "all";
my $sRunLog        = "rnaseq-submitted-jobs.txt";
my $sRunDir        = "./";
my $nRunCount      = 200;
my $sLsfAllocation = "acc_$ENV{USER}b";
my $nLsfMemory     = 50;
my $nLsfCores      = 12;
my $sLsfQueue      = "premium";
my $nLsfTime       = 24;
my $sLsfArch       = "bode,mothra";
GetOptions("help!"        => \$sHelp,
           "input:s"      => \$sInput,
           "genome:s"     => \$sGenome,
           "strand:s"     => \$sStrand,
           "task:s"       => \$sTask,
           "runcount:s"   => \$nRunCount,
           "log:s"        => \$sRunLog,
           "folder:s"     => \$sRunDir,
           "Project:s"    => \$sLsfAllocation,
           "queue:s"      => \$sLsfQueue,
           "memory:s"     => \$nLsfMemory,
           "cpu:s"        => \$nLsfCores,
           "architecture" => \$sLsfArch,
           "walltime:s"   => \$nLsfTime);

# PRINT HELP
$sHelp = 1 unless($sInput and $sGenome);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
   
   Script to run large numbers of RNA-seq pipeline mapping jobs at the same time.
   Handles different arguments per job as well as limiting the number of jobs
   queued at any given time.
   
   Arguments:
   -i --input <string>
     Tab-delimited input file
   -g --genome <string>
     Reference genome for mapping.
   -r --runcount <integer>
     Number of RNA-Seq mapping jobs to have in queue at any given time. Default: $nRunCount

   Optional:
   -s --strand <string>
     Global strand configuration, one of opposite, same or none. Default: $sStrand
     The global setting will be overridden by any strand definition in the input file.
   -t --task <string>
     Pipeline task to run. Default: $sTask
   -l --log <string>
     Name of runlog file with jobs that were already submitted to the queue. Default: $sRunLog 
   -f --folder <string>
     Name of folder containing subdirectories with fastq files.
     The fastq files in the folders need to have the same name
     as their parent folders, with .fastq.gz appended for single-
     end reads and _1.fastq.gz and _2.fastq.gz for paired-end reads. 
     Default: $sRunDir
   -P --project <string>
     Minerva allocation. Default: $sLsfAllocation
   -q --queue <string>
     Minerva LSF queue. Default: $sLsfQueue
   -a --architecture <string>
     Node architecture type (bode,mothra,amd). Default: $sLsfArch
   -m --memory <numeric>
     Amount of memory for each job in Gb. Default: $nLsfMemory
   -c --cpu <numeric>
     Number of CPU cores requested for each job. Default: $nLsfCores
   -w --walltime <numeric>
     Job walltime. Default: $nLsfTime

   -help
     This help message
   
HELP
}


##########
## MAIN ##
##########

# Check prerequisite binaries
die "Error: can't find bjobs in path\n" unless (can_run("bjobs"));
die "Error: can't find submitjob in path\n" unless (can_run("submitjob"));
die "Error: can't find run-rnaseq-pipeline in path\n" unless (can_run("run-rnaseq-pipeline"));

# Check arguments
$sStrand    = lc($sStrand);
$sLsfQueue  = lc($sLsfQueue);
my $rAllocs = get_allocations();
die "Error: cpu count must be numeric\n" unless ($nLsfCores =~ /^\d+$/);
die "Error: memory must be numeric\n" unless ($nLsfMemory =~ /^\d+$/);
die "Error: number of jobs must be numeric\n" unless ($nRunCount =~ /^\d+$/);
die "Error: strand must be either opposite, same or none\n" unless ($sStrand =~ /^(opposite|same|none)$/);
die "Error: allocation '$sLsfAllocation' does not exist\n" unless (exists $rAllocs->{$sLsfAllocation});
die "Error: queue '$sLsfQueue' does not exist\n" unless ($sLsfQueue =~ /^(alloc|expressalloc|low|premium|private)$/);

# Read entries from log file if it exists
my %hRunLog = ();
if (-e $sRunLog){
   open RUNLOG, $sRunLog or die "Error: can't open '$sRunLog': $!\n";
   while (<RUNLOG>){
      next if (/^\s*$/);
      next if (/^ *#/);
      s/[\n\r]+$//;
      $hRunLog{$_}++;
   }
   close RUNLOG;
}

# Read input file and mark entries that were already submitted before
my %hInput = ();
open INPUT, $sInput or die "Error: can't open '$sInput': $!\n";
my $sHeader = <INPUT>;
$sHeader =~ s/[\n\r]+$//;
my @asHeader = split /\t/, $sHeader;
die "Error: missing header in input file\n" unless ( lc(shift @asHeader) =~ "^file");
while (<INPUT>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sSample, @asFlags) = split /\t/, $_, -1;
   die "Error: duplicate entries found in input file\n" if (exists $hInput{$sSample});
   $hInput{$sSample}{submitted} = exists $hRunLog{$sSample} ? 1 : 0;  
   if (@asHeader){
      for ( my $i=0 ; $i<@asFlags ; $i++){
         $hInput{$sSample}{$asHeader[$i]} = $asFlags[$i];
      }
   }
}
close INPUT;


# Start wait loop for job submissions
my $cwd = cwd();
my $flLoop = 1;
while($flLoop){
   my $nQueuedJobs        = get_jobs();
   my $nSlots             = $nRunCount - $nQueuedJobs;
   my $nRemainingJobCount = 0;
   
   # Open log file for writing (append mode)
   open OUTLOG, ">>$sRunLog" or die "Error: can't append to run log: $!\n";
   warn "$nSlots job slots available\n";
   foreach my $sSample (keys %hInput){
      next if $hInput{$sSample}{submitted};
      if ($nSlots){
         if (-d "$sRunDir/$sSample"){
            chdir "$sRunDir/$sSample";
            my $sJobReads = "";
            
            # Check for paired or single-end read files
            if (-e "$sSample.fastq.gz"){
               $sJobReads = "single";
            }
            elsif ( (-e "${sSample}_1.fastq.gz") and (-e "${sSample}_2.fastq.gz") ){
               $sJobReads = "paired";
            }
            else{
               $sJobReads = "";
            }
            
            if ($sJobReads){
               my $sJobStrand = exists $hInput{$sSample}{strand} ? $hInput{$sSample}{strand} : $sStrand;
               my $sJobTask   = exists $hInput{$sSample}{task}   ? $hInput{$sSample}{task}   : $sTask;
               my $sJobGenome = exists $hInput{$sSample}{genome} ? $hInput{$sSample}{genome} : $sGenome;
               my $sJobQueue  = $sLsfQueue eq "private" ? "-q private -g bashia02" : "-q $sLsfQueue -P $sLsfAllocation";
               my $sJobArch   = $sLsfArch ? "-a $sLsfArch" : "";
               
               system("submitjob $nLsfTime -c $nLsfCores -m $nLsfMemory $sJobQueue $sJobArch run-rnaseq-pipeline $sJobReads $sJobGenome in=$sSample strand=$sJobStrand $sJobTask");
               print OUTLOG "$sSample\n";
               $hInput{$sSample}{submitted} = 1;
               $nSlots--;
            }
            else{
               warn("Skipping sample '$sSample' because there are no single or paired-end fastq files in the run folder\n");
            }
            chdir $cwd;
         }
         else{
            warn("Skipping sample '$sSample' because the sample dir does not exist\n");
         }
      }
      else{
         $nRemainingJobCount++;
      }
   }
   close OUTLOG;
   
   # Write out progress and see if there's anything left to be submitted
   print "$nRemainingJobCount jobs remaining to be submitted\n";
   if ($nRemainingJobCount){
      sleep 1800;
   }
   else{
      $flLoop = 0;
   }
}


#################
## SUBROUTINES ##
#################

# get_allocations
#
# Get hash of allocations
sub get_allocations {
   my %hAllocations;
   open BALANCE, "mybalance |" or die "Error: can't get cluster allocation listing: $!\n";
   while (<BALANCE>){
      $hAllocations{$1}++ if (/^(acc_\S+)\s+/);
   }
   close BALANCE;
   return \%hAllocations;
}

# get_jobs
#
# Get RNA-seq jobs in queue
sub get_jobs {
   open JOBS, "bjobs 2> /dev/null | grep run-rnaseq-p | wc -l |";
   my $nJobCount = <JOBS>;
   $nJobCount =~ s/[\n\r]+$//;
   close JOBS;
   return $nJobCount;
}

