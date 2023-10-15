#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use File::Temp qw(tempfile tempdir);

# GLOBALS
$ENV{TMPDIR}      ||= "/tmp";     # location for tmp file storage
$ENV{SORT}        ||= "sort";     # Unix sort binary
$ENV{SORT_BUFFER} ||= "2G";       # Maximum size of the sort buffer
$SIG{'INT'}=$SIG{'HUP'}=$SIG{'ABRT'}=$SIG{'QUIT'}=$SIG{'TRAP'}=$SIG{'STOP'}=\&INTERRUPT;

# ARGUMENTS
my $sHelp            = 0;
GetOptions("help!"   => \$sHelp);
my $sReadInput = shift @ARGV;           


# PRINT HELP
my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
$sHelp = 1 unless($sReadInput);
if ($sHelp) {
    die <<HELP

    $sScriptName <ordered-reads-file>

    Required:
    reads-file <string>
      Sequence read positions in basic bed format (chr, start, end, uniqueID)

    Options:
    -help
      This help message
      
HELP
}

# Check input
die "Error: file '$sReadInput' does not exist\n" unless (-e $sReadInput);

# Make sure the read input is really sorted properly
my $sResult = `$ENV{SORT} -s -c -t '\t' -k1,1 -k2,2n $sReadInput 2>&1`;
if ($sResult){
   chomp $sResult;
   warn "It appears that the reads file is not sorted:\n  ${sResult}\nstarting sort now...\n" if ($sResult);
   $sReadInput = sort_readinput($sReadInput);
}

# Open the reads file, find sets of overlapping reads and store the start and end positions
# Once we hit another chromosome or a non-overlapping area, process the overlapping area for bedgraph output.
my %hReads;
my ($sLastChr, $nLastEnd) = ('',0);
open IN, $sReadInput or die "Error: can't open '$sReadInput': $!\n";
print "track type=bed\n";
while (<IN>){
   next if /^\s*$/;
   next if /^\s*#/;
   s/[\n\r]$//g;
   my ($sChr, $nStart, $nEnd, $sID) = split /\t/;
   die "Error: read start position must be an integer on line $.\n" unless ($nStart =~ /^\d+$/);
   die "Error: read end position must be an integer on line $.\n"   unless ($nEnd =~ /^\d+$/);
   if ( ($sChr ne $sLastChr) or ($nStart > $nLastEnd) ){
      print_wigsegment($sLastChr, \%hReads) if (keys(%hReads));
      $sLastChr   = $sChr;
      $nLastEnd   = $nEnd;
      %hReads     = ();
      $hReads{$nStart}{'start'}++;
      $hReads{$nEnd}{'end'}++;
   }
   else{
      $nLastEnd = $nEnd if ($nEnd > $nLastEnd);
      $hReads{$nStart}{'start'}++;
      $hReads{$nEnd}{'end'}++;
   }
}
close IN;
# Still need to print the final line(s) of the file here
print_wigsegment($sLastChr, \%hReads) if (keys(%hReads));


#################
## SUBROUTINES ##
#################

# print_wigsegment
#
# Prints a wig segment of overlapping reads to STDOUT
sub print_wigsegment{
   my ($sChr, $rhReads) = @_;
   
   # Check the segment parameters before processing
   my @anPositions = sort {$a <=> $b} keys(%$rhReads);
   my $nSegStart = $anPositions[0];
   my $nSegEnd   = $anPositions[$#anPositions];
   die "Error: end position found at end of bed segment '$nSegStart-$nSegEnd\n"    if(exists($rhReads->{$nSegStart}{'end'}));
   die "Error: start position mission for bed segment '$nSegStart-$nSegEnd'\n"     unless(exists($rhReads->{$nSegStart}{'start'}));
   die "Error: start position found at end of bed segment '$nSegStart-$nSegEnd'\n" if(exists($rhReads->{$nSegEnd}{'start'}));
   die "Error: end position mission for bed segment '$nSegStart-$nSegEnd'\n"       unless(exists($rhReads->{$nSegEnd}{'end'}));
   
   # Process the hash, increment count for each start position
   # Print a bedgraph line when there is a change in position
   # If the last position was an end position, print the count and then
   # decrease the count by the number of end positions (not start!)
   my $nReadCount;
   my $nLastPosition;
   foreach my $nPosition (@anPositions){
      if ($nLastPosition){
         print join("\t", $sChr, $nLastPosition-1, $nPosition-1, $nReadCount), "\n";
         $nLastPosition = $nPosition;
         $nReadCount   += $rhReads->{$nPosition}{'start'} if (exists($rhReads->{$nPosition}{'start'}));
         $nReadCount   -= $rhReads->{$nPosition}{'end'}   if (exists($rhReads->{$nPosition}{'end'}));
      }
      else{
         $nLastPosition = $nPosition;
         $nReadCount   += $rhReads->{$nPosition}{'start'};
      }
   }
   die "Error: read count should be zero after bed segment processing, found '$nReadCount'\n" if ($nReadCount);
}


# sort_readinput
#
# Prepares a temporary file with sorted reads and returns the filename 
sub sort_readinput{
   my $sInput = shift @_;
   my ($fhTmp, $sTmp) = tempfile('sorted-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1);
   open READS, "$ENV{SORT} -S $ENV{SORT_BUFFER} -s -t '\t' -k1,1 -k2,2n $sInput |" or die "Error sorting reads file: $!\n";
   while (<READS>){
      print $fhTmp $_;
   }
   close READS;
   close $fhTmp;   
   return $sTmp;
}


# INTERRUPT
#
# Interrupt routine, make sure we exit gracefully for tmp file cleanup
sub INTERRUPT{
   exit(1); # This will call END
}
