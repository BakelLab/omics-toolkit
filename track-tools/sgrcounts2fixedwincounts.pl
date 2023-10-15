#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;

# GLOBALS
$ENV{SORT}        ||= "sort";     # Unix sort binary
$ENV{SORT_BUFFER} ||= "2G";       # Maximum size of the sort buffer

# ARGUMENTS
my $sHelp             = 0;
my $nWindowSize       = 100;
my $nMapPositions     = 0;
my $nSuppressZeroWins = 0;
GetOptions("help!"          => \$sHelp,
           "window_size:i"  => \$nWindowSize,
           "position!"      => \$nMapPositions,
           "suppress!"      => \$nSuppressZeroWins);
my $sSgrFile = shift @ARGV;           


# PRINT HELP
my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
$sHelp = 1 unless($sSgrFile);
if ($sHelp) {
    die <<HELP

    $sScriptName [ -w -s -p ] <sgr-file>

    Convert an sgr file with count data at mapped positions in the genome
    to a file with counts in non-overlapping windows of a specified width.
    The input file must be sorted by chromosome and position.

    Options:
    -w <integer>
      The window size in bp
      default: $nWindowSize
    -s
      Supress printing windows without any read coverage
    -p
      Output middle position of sequence window instead of start and end
    -help
      This help message
      
HELP
}

# Check input
unless (-e $sSgrFile)      {die "File '$sSgrFile' does not exist\n"}
unless ($nWindowSize >  0 )  {die "Invalid option for -window_size: '$nWindowSize'\n";}

# Make sure the read input is really sorted properly
my $sResult = `$ENV{SORT} -s -c -t '\t' -k1,1 -k2,2n $sSgrFile 2>&1`;
die "Error: reads file is not sorted by chromosome and start position\n$sResult\n" if ($sResult);

# Now sort the filter file and process to get coverage and read count
&process_filtered_reads($sSgrFile, $nWindowSize, $nMapPositions, $nSuppressZeroWins);


#################
## SUBROUTINES ##
#################

# Get read counts in fixed windows across the genome
sub process_filtered_reads{
   my ($sReadFile, $nWindowSize, $nMapPositions, $nSuppressZeroWins) = @_;
   
   # Process the sorted reads file
   open READS, $sReadFile or die "Can't open '$sReadFile': $!\n";
   my %hWindows;
   my %hLastPrintWindow;
   while (<READS>){
      my ($sChr, $nPos, $nCount) = (split /\t/)[0..2];
      $hLastPrintWindow{$sChr} = -1 unless(exists($hLastPrintWindow{$sChr}));
      my $nReadWindow  = int($nPos/$nWindowSize);
      unless (exists $hWindows{$sChr}{$nReadWindow}){
         &print_windows(\%hWindows, \%hLastPrintWindow, $nWindowSize, $nMapPositions, $nSuppressZeroWins);
         %hWindows = ();
      }
      $hWindows{$sChr}{$nReadWindow} += $nCount; 
   }
   &print_windows(\%hWindows, \%hLastPrintWindow, $nWindowSize, $nMapPositions, $nSuppressZeroWins);
   close READS;
}

# Print the windows that are currently in the window hash
sub print_windows{
   my ($rhWindows, $rhLastPrintWindow, $nWindowSize, $nMapPositions, $nSuppressZeroWins) = @_;
   foreach my $sChr (sort keys(%$rhWindows)){
      foreach my $nWindow (sort {$a <=> $b} keys %{$rhWindows->{$sChr}}){
         # Print previous windows that didn't have any coverage
         unless ($nSuppressZeroWins){
            for (my $i=$rhLastPrintWindow->{$sChr}+1 ; $i<$nWindow ; $i++){
               if ($nMapPositions){
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
         my $nCount = $rhWindows->{$sChr}{$nWindow};
         if ($nMapPositions){
            my $nPosition = ($nWindow * $nWindowSize) + int($nWindowSize/2);
            print join("\t", $sChr, $nPosition, $nCount), "\n";
         }
         else{
            my $nStart = ($nWindow * $nWindowSize) + 1;
            my $nEnd   = ($nWindow * $nWindowSize) + $nWindowSize;
            print join("\t", $sChr, $nStart, $nEnd, $nCount), "\n";
         }
      }
   }
}
