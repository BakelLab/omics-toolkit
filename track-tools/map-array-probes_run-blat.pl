#!/usr/bin/env perl

# map-array-probes_run-blat
# Maps microarray probes to a genome using BLAT

# MODULES
use strict;
use Getopt::Long;
use File::Basename;

# ENVIRONMENT
$ENV{'BLAT_BINARY'}   = 'blat' unless $ENV{'BLAT_BINARY'};
$ENV{'QUEUE_WRAPPER'} = 'submitjob' unless $ENV{'QUEUE_WRAPPER'};

# GET PARAMETERS
my $sHelp     = 0;
my $sInput    = '';
my $sOutdir   = '';
my $nCpuTime  = 48;
GetOptions("help!"     => \$sHelp,
           "input:s"   => \$sInput,
           "output:s"  => \$sOutdir,
           "time:i"    => \$nCpuTime);

# PRINT HELP
$sHelp = 1 unless($sInput and $sOutdir and @ARGV and $nCpuTime);
if ($sHelp) {
    die <<HELP
    map-array-probes_run-blat -i <probe sequence file> -o <output dir> genome-seq-1.fa ... genome-seq-n.fa
    
    Takes a tab-delimited file containing probe identifiers and sequence
    and maps it against a set of fasta formatted genome files using blat.
    Each blat search of a genome file will be submitted to the cluster as a
    separate job.
    To parse the blat output, run the map-array-probes_parse-blat.pl script.
    
    arguments:
    -i Input file
        Tab delimited list of probe ID and sequence
         probeID <tab> sequence
    -o Output dir
        The output dir that will hold all blat output files
    -t Expected maximum run time
        Maximum expected time for a blat search against the largest chromosome in
        the list of fasta files. This time will be used when submitting to the 
        torque queue.
        Default: $nCpuTime hours
    -h Help
HELP
}


######################################
#               START                #
######################################

# Create output dir if it doesn't exist yet
unless (-e $sOutdir) { mkdir($sOutdir) or die "Could not create mapping output dir: $!\n"};

# First reformat the tab-delimited input file into a temporary fasta file
print "Preparing sequence file for blat\n";
my $sInputFasta = &create_fasta_file($sInput, $sOutdir);

# Now start a blat for each genome file
print "Submitting jobs to queue\n";
foreach my $sGenomeFile (@ARGV){
   if (-e $sGenomeFile){
      print "Starting blat job for $sGenomeFile\n";
      &start_blat_job($sInputFasta, $sGenomeFile, $sOutdir, $nCpuTime);
   }
   else{
      print STDERR "$sGenomeFile does not exist, skipping!\n";
   }
}


######################################
#            SUBROUTINES             #
######################################

# unique_name()
#
# Returns a unique string of 15 characters
sub unique_name {
    my @asCh   = ( "A" .. "Z", 0 .. 9 );
    my $sName  = join( "", @asCh[ map { rand @asCh } ( 1 .. 5 ) ], time() );
    return $sName;
}


# create_fasta_file()
#
# Converts the tab-delimited input file in a fasta formatted input file
sub create_fasta_file {
   my ($sInput, $sOutdir) = @_;
   my $sTmpFasta = join('', 'temp-', &unique_name(), '.fa');
   
   open TMPFASTA, ">$sOutdir/$sTmpFasta" or die "Can't create temporary fasta input file: $!\n";
   open INPUT, $sInput or die "Can't read from probe sequence file: $!\n";
   while (<INPUT>){
      next if (/^\s*$/);
      next if (/^\s*\t*#/);
      next if (/^ProbeID/);
      s/[\n\r]//g;
      my ($sID, $sSeq) = split /\t/;
      print TMPFASTA join('', '>', $sID, "\n", $sSeq, "\n");
   }
   close TMPFASTA;
   close INPUT;

   return($sTmpFasta);
}


# start_blat_job()
#
# Submits a blat job to the cluster
sub start_blat_job {
   my ($sInputFasta, $sGenomeFile, $sOutdir, $nCpuTime) = @_;
   
   # Create PBS output dir if it doesn't exist yet
   my $sPBSout = "$ENV{HOME}/pbs-output";
   unless (-e $sPBSout) { mkdir($sPBSout) or die "Could not create pbs output dir: $!\n"};
   
   # Get genome file basename
   my $sGenomeBase = basename($sGenomeFile);
   
   # Submit to queue
   `$ENV{QUEUE_WRAPPER} $nCpuTime $ENV{'BLAT_BINARY'} $sGenomeFile $sOutdir/$sInputFasta $sOutdir/$sGenomeBase.psl -tileSize=7 -minScore=18 -noHead`
}
