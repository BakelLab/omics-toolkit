#!/usr/bin/perl

# 15.12.2010 13:14:46 EST
# Harm van Bakel <hvbakel@gmail.com>

# GENERAL
my $sRunJob    = '/home/vanbakel/bin/run-overlapper-jobs.mk';   # Location of the runjob makefile
my $sSubmitjob = 'submitjob';                                   # Location of the submitjob script

# MODULES
use strict;
use warnings;
use Getopt::Long;
use Cwd;

# GET PARAMETERS
my $sHelp        = 0;
my $nBstart      = 1;
my $nBend        = 5000;
my $nConcurrency = 120;
my $nThreads     = 1;
GetOptions("help!"         => \$sHelp,
           "bstart:n"      => \$nBstart,
           "bend:n"        => \$nBend,
           "concurrency:n" => \$nConcurrency,
           "threads:n"     => \$nThreads);

# PRINT HELP
$sHelp = 1 unless($nBstart and $nBend and $nConcurrency and $nThreads);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
    
    -bstart <integer>
      The overlap batch ID to start at. Default: $nBstart
    -bend <integer>
      The overlap batch ID to end at. Default: $nBend
    -concurrency <integer>
      The maximum number of jobs to run at any given time
      default: $nConcurrency
    -threads <integer>
      The number of overlapper jobs to run per node
      default: $nThreads
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check arguments
my @asErrors;
push @asErrors, "Argument 'bstart' must be a positive integer"      unless ($nBstart =~ /^[1-9]\d*$/);
push @asErrors, "Argument 'bend' must be a positive integer"        unless ($nBend =~ /^[1-9]\d*$/);
push @asErrors, "Argument 'concurrency' must be a positive integer" unless ($nConcurrency =~ /^[1-9]\d*$/);
push @asErrors, "Argument 'threads' must be a positive integer"     unless ($nThreads =~ /^[1-9]\d*$/);

# Make sure the script is run in the proper folder and that all prerequisites are available
push @asErrors, "The 'runjob.mk' makefile doesn't exist at the indicated location ($sRunJob)" unless (-e $sRunJob);
push @asErrors, "This script must be executed in an 'overlap' folder"             unless (cwd() =~ /overlap/);
push @asErrors, "The 'ovlbat' file was not found in the current folder"     unless (-e 'ovlbat');
push @asErrors, "The 'ovljob' file was not found in the current folder"     unless (-e 'ovljob');
push @asErrors, "The 'ovlopt' file was not found in the current folder"     unless (-e 'ovlopt');
push @asErrors, "The 'overlap.sh' file was not found in the current folder" unless (-e 'overlap.sh');

# Output if any errors were found
if (@asErrors){
   my $sErrorMsg = "The following errors were found:\n";
   foreach my $sError (@asErrors){$sErrorMsg .= "   $sError\n"}
   die $sErrorMsg;
}

# Adjust bend to the max job count
my $nMaxJobCount = get_job_count('ovljob');
if ($nBend > $nMaxJobCount){
   warn("Adjusting 'bend' from $nBend to $nMaxJobCount (the maximum number of overlap jobs)\n");
   $nBend = $nMaxJobCount;
}

# Now we can start the job submission fun
my %hDepend;
my $nDependCount = 1;
my $i = $nBstart;
while ($i <= $nBend){
   my @asBatch;
   my $sJobID    = "";
   $nDependCount = 1 if ($nDependCount > $nConcurrency);
   
   
   for (my $j=0 ; $j<$nThreads ; $j++){push @asBatch, $i++};
   my $sBatch = join(' ', @asBatch);
   if (exists $hDepend{$nDependCount}){
      $sJobID = `$sSubmitjob 8 -c 8 -m 11 -W depend=afterany:$hDepend{$nDependCount} $sRunJob -j $nThreads batch='$sBatch' 2>&1 1>/dev/null`;
   }
   else{
      $sJobID = `$sSubmitjob 8 -c 8 -m 11 $sRunJob -j $nThreads batch='$sBatch' 2>&1 1>/dev/null`;
   }
   $sJobID =~ s/[\n\r]+$//;
   die "Error: job submission failed (qsub returned '$sJobID')\n" unless ($sJobID =~ /^\d+\..*$/);
   $hDepend{$nDependCount} = $sJobID;
   $nDependCount++;
   print "$sJobID\n";
}

#################
## SUBROUTINES ##
#################

# get_job_count
#
# Count the number of overlap jobs
sub get_job_count {
   my $sFile = shift @_;
   my $nMaxJobCount = 0;
   open COUNT, $sFile or die "Error: can't open '$sFile': $!\n";
   while (<COUNT>){
      next if (/^\s*$/);
      next if (/^ *#/);
      s/[\n\r]+$//;
      s/^0+//;
      $nMaxJobCount = $_ if ($_ > $nMaxJobCount);
   }
   close COUNT;
   return $nMaxJobCount;
}
