#!/usr/bin/env perl

# 11.06.2013 12:07:46 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use String::Approx 'amatch';
use IO::File;
use IO::Zlib qw(:gzip_external 1);

# GET PARAMETERS
my $sHelp         = 0;
my $sFastq        = "";
my $sMatchFile    = "";
my $nEditDist     = 0;
my $flRevComp     = 0;
my $flInvertMatch = 0;
GetOptions ("fastq=s"         => \$sFastq,
            "match=s"         => \$sMatchFile,
            "editdist:n"      => \$nEditDist,
            "help!"           => \$sHelp,
            "revcomp!"        => \$flRevComp,
            "v|invert-match!" => \$flInvertMatch) || die "Bad option";

# PRINT HELP
$sHelp = 1 unless($sFastq and ($sMatchFile or @ARGV));
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -f <fastq file> [-m <match-file>] [<seq1> ... <seqN>]
    
   Filter reads in a fastq file based on sequence matches
   
   Arguments:
    -f --fastq <string>
      Name of the file with fastq data
    -m --match <string>
      Name of the file with sequences to match
    -e --editdist <integer>
      Edit distance (sum of the number of allowed insertions, deletions,
      and substitutions). Default: $nEditDist
    -r --revcomp
      Match reverse-complement as well
    -v --invert-match
      Invert the sense of matching, to select non-matching reads.
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Read the sequences to filter
my %hMatchSeqs;
if ($sMatchFile){
   open IN, $sMatchFile or die "Error\n";
   while (<IN>){
      next if (/^\s*$/);
      next if (/^ *#/);
      my $sMatchSeq = $_;
      $sMatchSeq =~ s/[\n\r]+$//;
      $hMatchSeqs{$sMatchSeq} = 0;
   }
   close IN;
}

# Read sequences from the command line
foreach my $sMatchSeq (@ARGV){
   $hMatchSeqs{$sMatchSeq} ||= 0;
}

# Open the input files
my ($fhFastq);
if($sFastq =~ /\.gz$/) {
   $fhFastq = new IO::Zlib($sFastq, "rb") or die "Error: can't open '$sFastq': $!\n";
}
else{
   $fhFastq = new IO::File($sFastq, "r") or die "Error: can't open '$sFastq': $!\n";
}

# Process the reads
my $nExtracted = 0;
while(<$fhFastq>) {
   my $name = $_;
   my $seq  = <$fhFastq>;
   my $tmp  = <$fhFastq>;
   my $qual = <$fhFastq>;
   if ($qual){
      die "Error: read name '$name' in file '$sFastq' does not conform to the fastq standard\n" unless ($name =~ /^\@/);
      my @seq = ($seq);
      push @seq, reverse_complement($seq) if ($flRevComp);
      my $flExtract = 0;
      foreach my $sMatchSeq (keys %hMatchSeqs){
         if (amatch($sMatchSeq, ["i $nEditDist"], (@seq))){
            $flExtract = 1;
            last;
         }
      }
      
      $flExtract = 1 - $flExtract if ($flInvertMatch);
      if ($flExtract){
         print "$name$seq$tmp$qual";
         $nExtracted++;
      }
   }
   else{
      die "Error: Missing reads at end of file\n";
   }
}
close($fhFastq);

# Report the number of sequences extracted
my $sMessage = $nExtracted == 1 ? "Extracted $nExtracted read\n" : "Extracted $nExtracted reads\n";
warn ($sMessage);


###############
# SUBROUTINES #
###############

# reverse_complement
#
# Returns the reverse-complement of the supplied sequence
sub reverse_complement{
   my $seq = shift(@_);
   my $rev = reverse $seq;
   $rev =~ tr/ACGTacgt/TGCAtgca/;
   return $rev;
}
