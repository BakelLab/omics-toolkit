#!/usr/bin/env perl

use strict;
use Getopt::Long;

# GET ARGUMENTS
my $sHelp       = 0;
my $nOverlapCut = 0.5;
my $nBuffer     = 350;
GetOptions("help!"      => \$sHelp,
	   "overlap:s"  => \$nOverlapCut,
	   "buffer:i"   => \$nBuffer);

# PRINT HELP
$sHelp = 1 unless($nOverlapCut and @ARGV>1);
if ($sHelp) {
    die <<HELP
    
    peakmatcher.pl [-o overlap] peakfinderA.csv peakfinderB.csv
    Find overlap between peaks in the two files supplied
    
    arguments:
    -o Overlap
        Give the minimum fraction of overlap for a peakmatch
        Default: 0.5
    -b Buffer size
        Give a buffer size within which the peaks are still considered
        significantly overlapping
        Default: 350
    -h Help
    
HELP
}


######################################
#               START                #
######################################

my $sPeakfinderA = shift @ARGV;
my $sPeakfinderB = shift @ARGV;

# Hash the first input file
my %hPeakfinderA;
my $nPeakfinderA_total; 
open INPUT, $sPeakfinderA or die "Can't open $sPeakfinderA: $!\n";
while (<INPUT>){
    s/[\n\r]//g;
    next if (/^\s*$/);
    next if (/^\s*#/);
    my ($sFeatures, $sCall, $nCallPerc, $sChr, $nStart, $nStop, $nLength, $nOrfLength, $nIGlength, $nNoOligos, $nAvgRatio, $nMaxRatio, $nMinP, $nIgbUrl) = split /\t/;
    push @{$hPeakfinderA{$sChr}}, [$nStart, $nStop];
    $nPeakfinderA_total++;
}
close INPUT;

# Now match with the second input file
my $nPeakfinderB_overlap;
my $nPeakfinderB_unique;
my @aaOverlapList;
my @aaPeakBList;
open INPUT, $sPeakfinderB or die "Can't open $sPeakfinderA: $!\n";
while (<INPUT>){
    s/[\n\r]//g;
    next if (/^\s*$/);
    next if (/^\s*#/);
    my ($sFeatures, $sCall, $nCallPerc, $sChr, $nStart, $nStop, $nLength, $nOrfLength, $nIGlength, $nNoOligos, $nAvgRatio, $nMaxRatio, $nMinP, $nIgbUrl) = split /\t/;
    my $nOverlap = &get_peak_overlap($nStart, $nStop, $hPeakfinderA{$sChr}, $nBuffer) if (exists($hPeakfinderA{$sChr}));
    if ($nOverlap){ 
	$nPeakfinderB_overlap++;
	push @aaOverlapList, [$sFeatures, $sChr, $nStart, $nStop, $nMaxRatio];
    }
    else { 
	$nPeakfinderB_unique++;
	push @aaPeakBList, [$sFeatures, $sChr, $nStart, $nStop, $nMaxRatio];
    }
}
close INPUT;

# Print overlap numbers
print "Unique to $sPeakfinderA\t", $nPeakfinderA_total-$nPeakfinderB_overlap, "\n";
print "Unique to $sPeakfinderB\t", $nPeakfinderB_unique, "\n";
print "Peak overlap $sPeakfinderB\t", $nPeakfinderB_overlap, "\n";

# Print overlap lists for the B group
print "\n\nOverlapping peaks in $sPeakfinderB:\n";
foreach my $rOverlap (@aaOverlapList){
    print join ("\t", @$rOverlap), "\n";
}

print "\n\nUnique peaks in $sPeakfinderB:\n";
foreach my $rPeakB (@aaPeakBList){
    print join ("\t", @$rPeakB), "\n";
}


######################################
#           SUBROUTINES              #
######################################


sub get_peak_overlap {
    my ($nBstart, $nBstop, $rApeaks, $nBuffer) = @_;
    my $nReturn = 0;
    $nBuffer = $nBuffer/2;
    $nBstart = $nBstart - $nBuffer;
    $nBstop  = $nBstop  + $nBuffer;
    
    foreach my $rApeak (@$rApeaks){
	my ($nAstart, $nAstop) = @$rApeak;
	$nAstart = $nAstart - $nBuffer;
	$nAstop  = $nAstop  + $nBuffer;
	my $nMinsize  = $nAstop-$nAstart+1;
	my $nMaxStart = $nAstart;
	my $nMinStop  = $nAstop;
	$nMinsize     = $nBstop-$nBstart+1 if ( ($nBstop-$nBstart+1) <$nMinsize);
	$nMaxStart    = $nBstart if ($nBstart>$nMaxStart);
	$nMinStop     = $nBstop  if ($nBstop<$nMinStop);
	my $nOverlap  = $nMinStop-$nMaxStart+1;
	
	if ($nOverlap>0){
	    $nReturn = 1 if ( ($nOverlap/$nMinsize) > $nOverlapCut);
	}
    }
    
    return ($nReturn);
}




