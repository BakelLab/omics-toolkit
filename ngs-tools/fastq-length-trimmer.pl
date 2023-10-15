#!/usr/bin/env perl

# Modified from a bowtie script written by Ben Langmead

# MODULES
use strict;
use warnings;
use Getopt::Long;
use IO::File;
use IO::Zlib qw(:gzip_external 1);

# GET PARAMETERS
my $sHelp         = 0;
my $nTrimLeft     = 0;
my $nKeepRight    = 0;
my $nMaxLength    = 0;
my $flAppendRange = 0;
GetOptions ("maxlength:n"  => \$nMaxLength,
            "trimleft:n"   => \$nTrimLeft,
            "keepright:n"  => \$nKeepRight,
            "appendrange!" => \$flAppendRange) || die "Bad option";
my $sInput = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sInput);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [-tkm] <fastq-file>
    
   Trim reads in a fastq input file to a specified length by removing
   bases at the 5' and/or 3' ends. Reads uncompressed or gzip compressed files, 
   depending on whether a .gz extension is present in the filename
   
   Arguments:
    -t --trimleft <integer>
      Number of bases to trim from the 5' end of each read.
      default: no trimming
    -k --keepright <integer>
      Trim the 5' end of each read, keeping this number of 3' bases
      default: no trimming
    -m --maxlength <integer>
      Trim all reads at the 3' end to this maximum length (applied *after*
      5' end trimming).
      default: no trimming
    -a --appendrange <integer>
      Append the trimmed range to the read identifier
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########


# Open the input files
my ($fhInput);
if($sInput =~ /\.gz$/) {
   $fhInput = new IO::Zlib($sInput, "rb") or die "Error: can't open '$sInput': $!\n";
}
else{
   $fhInput = new IO::File($sInput, "r") or die "Error: can't open '$sInput': $!\n";
}

# Check arguments
die "Error: the -t and -k arguments cannot be used together\n" if ($nKeepRight and $nTrimLeft);

# Process the reads
my $nSkipped = 0;
while(<$fhInput>) {
   my $sName   = $_;
   my $sSeq    = <$fhInput>;
   my $sTmp    = <$fhInput>;
   my $sQual   = <$fhInput>;
   $sSeq  =~ s/[\n\r]+$//;
   $sQual =~ s/[\n\r]+$//;
   die "Error: mismatch between length of sequence and quality string on line $.\n" unless(length($sSeq) == length($sQual));
   my $nSeqLen = length($sSeq);
   my ($nLeftTrim, $nRightTrim) = (0,0);
   if ($nTrimLeft){
      my $nPreTrim = length($sSeq);
      $sSeq  = $nPreTrim > $nTrimLeft ? substr($sSeq,$nTrimLeft)  : "";
      $sQual = $nPreTrim > $nTrimLeft ? substr($sQual,$nTrimLeft) : "";
      $nLeftTrim = $nPreTrim - length($sSeq);
   }
   if ($nKeepRight){
      my $nPreTrim = length($sSeq);
      $sSeq  = $nPreTrim > $nKeepRight ? substr($sSeq,  -$nKeepRight) : $sSeq;
      $sQual = $nPreTrim > $nKeepRight ? substr($sQual, -$nKeepRight) : $sQual;
      $nLeftTrim = $nPreTrim - length($sSeq);
   }
   if ($nMaxLength){
      my $nPreTrim = length($sSeq);
      $sSeq  = substr($sSeq,0,$nMaxLength)  if ($nPreTrim > $nMaxLength);
      $sQual = substr($sQual,0,$nMaxLength) if ($nPreTrim > $nMaxLength);
      $nRightTrim = $nPreTrim - length($sSeq);
   }
   if ($flAppendRange){
      $sName =~ s/[\n\r]+$//;
      my $nStart = $nLeftTrim+1;
      my $nEnd   = $nSeqLen - $nRightTrim;
      if ($sName =~ /(^.*)(\/\d)$/){
         $sName = join('', $1, ":$nStart-$nEnd", $2) . "\n";
      }
      elsif ($sName =~ /(^.*) (.*)$/){
         $sName = join('' , $1, ":$nStart-$nEnd ", $2) . "\n";
      }
      else{
         $sName = join(':', $sName, "$nStart-$nEnd") . "\n";
      }
   }
   print $sName,$sSeq,"\n",$sTmp,$sQual,"\n";
}
close $fhInput;
