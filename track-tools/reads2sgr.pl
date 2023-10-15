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
my $nWindowSize      = 40;
my $flBedBasic       = 0;
my $flPrintZeroWins  = 0;
my $flCoverage       = 0;
my $flSingleCount    = 0;
GetOptions("help!"          => \$sHelp,
           "window_size:i"  => \$nWindowSize,
           "bedbasic!"      => \$flBedBasic,
           "zerowins!"      => \$flPrintZeroWins,
           "coverage!"      => \$flCoverage,
           "single!"        => \$flSingleCount);
my $sReadInput = shift @ARGV;           


# PRINT HELP
my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
$sHelp = 1 unless($sReadInput);
if ($sHelp) {
    die <<HELP

    $sScriptName [ -w -s -p ] <ordered-reads-file>

    Required:
    reads-file <string>
      Sequence read positions in basic bed format (chr, start, end, uniqueID)

    Options:
    -w <integer>
      The window size in bp
      default: $nWindowSize
    -z
      Print windows without any read coverage, default is to skip windows without reads
    -b
      Output window start and end (bed_basic) instead of just the middle position (sgr)
    -c
      Output base coverage instead of read count in each window
    -s
      Make sure each read is only counted once (window assignment based on read midpoint)
    -help
      This help message
      
HELP
}

# Check input
unless (-e $sReadInput)      {die "File '$sReadInput' does not exist\n"}
unless ($nWindowSize >  0 )  {die "Invalid option for -w: '$nWindowSize'\n";}

# Make sure the read input is really sorted properly
my $sResult = `$ENV{SORT} -s -c -t '\t' -k1,1 -k2,2n $sReadInput 2>&1`;
if ($sResult){
   chomp $sResult;
   warn "It appears that the reads file is not sorted:\n  ${sResult}\nstarting sort now...\n" if ($sResult);
   $sReadInput = sort_readinput($sReadInput);
}

# Now sort the filter file and process to get coverage and read count
&process_filtered_reads($sReadInput, $nWindowSize, $flBedBasic, $flPrintZeroWins, $flSingleCount);


#################
## SUBROUTINES ##
#################

# Round number to integer
sub round {
    my($nNumber) = shift;
    return int($nNumber + .5);
}

# Get read counts in fixed windows across the genome
sub process_filtered_reads{
   my ($sReadFile, $nWindowSize, $flBedBasic, $flPrintZeroWins, $flSingleCount) = @_;
   
   # Process the sorted reads file
   open READS, $sReadFile or die "Can't open '$sReadFile': $!\n";
   my %hWindows;
   my %hLastPrintWindow;
   while (<READS>){
      my ($sChr, $nStart, $nEnd) = (split /\t/)[0..2];
      ($nStart, $nEnd) = ($nEnd, $nStart) if ($nEnd < $nStart);
      my $nMid         = $nStart + round(($nEnd-$nStart+1)/2);
      $hLastPrintWindow{$sChr} = -1 unless(exists($hLastPrintWindow{$sChr}));
      my $nReadWindow       = int($nStart/$nWindowSize);
      my $nReadCenterWindow = int($nMid/$nWindowSize);
      unless (exists($hWindows{$sChr}{$nReadWindow})){
         &print_windows(\%hWindows, \%hLastPrintWindow, $nWindowSize, $flBedBasic, $flPrintZeroWins);
         %hWindows = ();             # clear the hash
      }
      while(my $nCoverage = &get_window_coverage($nStart, $nEnd, $nReadWindow, $nWindowSize)){
         $hWindows{$sChr}{$nReadWindow}{'coverage'} += $nCoverage;
         if ($flSingleCount){
				if ($nReadWindow == $nReadCenterWindow){  
				   $hWindows{$sChr}{$nReadWindow}{'count'}+= 1;
				}
				else{
				   $hWindows{$sChr}{$nReadWindow}{'count'}+= 0;
				}
			}
			else{
				$hWindows{$sChr}{$nReadWindow}{'count'}    += 1;
			}
         $nReadWindow++;
      }
   }
   &print_windows(\%hWindows, \%hLastPrintWindow, $nWindowSize, $flBedBasic, $flPrintZeroWins);
   close READS;
}


# Calculates the coverage of a read in a particular window
sub get_window_coverage{
   my ($nReadStart, $nReadEnd, $nReadWindow, $nWindowSize) = @_;
   my $nWindowStart  = ($nReadWindow * $nWindowSize) + 1;
   my $nWindowEnd    = $nWindowStart + $nWindowSize - 1;
   my $nOverlapStart = $nReadStart > $nWindowStart ? $nReadStart : $nWindowStart;
   my $nOverlapEnd   = $nReadEnd   < $nWindowEnd   ? $nReadEnd   : $nWindowEnd;
   my $nOverlap      = ($nOverlapEnd - $nOverlapStart + 1) <= 0 ? 0 : $nOverlapEnd - $nOverlapStart + 1;
   return ($nOverlap / $nWindowSize);
}


# Print the windows that are currently in the window hash
sub print_windows{
   my ($rhWindows, $rhLastPrintWindow, $nWindowSize, $flBedBasic, $flPrintZeroWins) = @_;
   foreach my $sChr (sort keys(%$rhWindows)){
      foreach my $nWindow (sort {$a <=> $b} keys %{$rhWindows->{$sChr}}){
         # Print previous windows that didn't have any coverage
         if ($flPrintZeroWins){
            for (my $i=$rhLastPrintWindow->{$sChr}+1 ; $i<$nWindow ; $i++){
               unless ($flBedBasic){
                  my $nPosition = ($i * $nWindowSize) + int($nWindowSize/2);
                  print join("\t", $sChr, $nPosition, 0), "\n";
               }
               else{
                  my $nStart = ($i * $nWindowSize) + 1;
                  my $nEnd   = ($i * $nWindowSize) + $nWindowSize;
                  print join("\t", $sChr, $nStart, $nEnd, 0), "\n";
               }
            }
         }
         $rhLastPrintWindow->{$sChr} = $nWindow;
      
         # Print the current window that has read coverage
         my ($nCount, $nCoverage) = ($rhWindows->{$sChr}{$nWindow}{'count'}, $rhWindows->{$sChr}{$nWindow}{'coverage'});
         my $nValue    = $flCoverage ? $nCoverage : $nCount;
         unless ($flBedBasic){
            my $nPosition = ($nWindow * $nWindowSize) + int($nWindowSize/2);
            print join("\t", $sChr, $nPosition, $nValue), "\n";
         }
         else{
            my $nStart = ($nWindow * $nWindowSize) + 1;
            my $nEnd   = ($nWindow * $nWindowSize) + $nWindowSize;
            print join("\t", $sChr, $nStart, $nEnd, $nValue), "\n";
         }
      }
   }
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
