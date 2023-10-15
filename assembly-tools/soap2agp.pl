#!/usr/bin/env perl

# 07.09.2010 20:00:48 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use FileHandle;

# GET PARAMETERS
my $sHelp         = 0;
my $nMinContig    = 200;
my $nMinGap       = 10;
GetOptions("help!"       => \$sHelp,
           "mincontig:n" => \$nMinContig,
           "mingap:n"    => \$nMinGap);
my $sInput  = shift @ARGV;
my $sOutput = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sInput and $sOutput);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [options] <file.fa> <output prefix>
   
   Convert a soap file to AGP/FSA files suitable for submission
   to genbank.
    
   Options:
    -mincontig <integer>
      Minimum allowed contig length. Contigs shorter than this
      cutoff will be marked as gaps and dropped from the assembly.
      Default: $nMinContig
    -mingap <integer>
      Minimum gap size. Smaller gaps will be increased to this length.
      default: $nMinGap
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Process the fasta file
my $sFastaHeader = '';
my $sFastaSeq    = '';
my $sLastSeq     = '';
open INPUT, "<$sInput" or die "Error: can't read '$sInput': $!\n";
my $fhOutAGP = new FileHandle;
my $fhOutFSA = new FileHandle;
$fhOutAGP->open(">$sOutput.agp") or die "Error: can't write to '$sOutput.agp': $!\n";
$fhOutFSA->open(">$sOutput.fsa") or die "Error: can't write to '$sOutput.fsa': $!\n";
while (<INPUT>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   if (/^>/ or eof){
      if (eof){
         die "Error: file ends in fasta header without sequence\n" if (/^>/);
         $sFastaSeq .= $_;
      }
      if ($sFastaHeader){
         my ($sID, @asRest) = split / +/, $sFastaHeader;
         $sID  =~ s/^>//;
         print_agp($sID, $sFastaSeq, $nMinContig, $fhOutAGP, $fhOutFSA);
      }
      $sFastaHeader = $_;
      $sFastaSeq    = "";
   }
   else{
      next if (/^\s*$/);
      next if (/^ *#/);
      $sFastaSeq .= $_ if ($sFastaHeader);
   }
}
close INPUT;
$fhOutAGP->close();
$fhOutFSA->close();


###############
# SUBROUTINES #
###############

# print_agp
#
# print agp fasta entry
sub print_agp {
   my ($sID, $sSeq, $nMinLen, $fhAGP, $fhFSA) = @_;
   
   # Only process sequences larger than the minimum length
   if (length($sSeq) >= $nMinLen){
   
      # Split sequence in segments based on 'N' characters and merge contigs
      # that do not meet the minimum size threshold into the neighboring gap.
      $nMinLen--;
      $sSeq =~ s/([nN]+)/:$1:/g;              # Delineate any segment with N's (i.e. gaps) with ':'
      $sSeq =~ s/^([^nN]{1,$nMinLen}):/$1/;   # Check first contig length
      $sSeq =~ s/:([^nN]{1,$nMinLen})$/$1/;   # Check last contig length
      $sSeq =~ s/:([^nN]{1,$nMinLen}):/$1/g;  # Check internal contig lengths
      my (@asSegments) = split /:/, $sSeq;    # Array of segments; gaps will contain at least one 'N'
      
      # Process each segment into an agp and corresponding fsa entry
      my $nContigCount = 0;
      my $nObjectCount = 0;
      my $nObjectStart = 0;
      my $nObjectEnd   = 0;
      for (my $i=0 ; $i<@asSegments ; $i++){
         next if ($asSegments[$i] =~ /[nN]/ and $nObjectCount == 0); # Make sure we don't start with a gap
         next if ($asSegments[$i] =~ /[nN]/ and $i == $#asSegments); # Make sure we don't end on a gap
         $asSegments[$i] = "N"x$nMinGap if ($asSegments[$i] =~ /[nN]/ and length($asSegments[$i])<$nMinGap); # Pad small gaps to minimum size
         $nObjectCount++;
         $nObjectStart = $nObjectEnd + 1;
         $nObjectEnd   = $nObjectEnd + length($asSegments[$i]);
         if ($asSegments[$i] =~ /[nN]/ ){   # we have a gap, print a gap line to the agp file
            print $fhAGP join("\t", $sID, $nObjectStart, $nObjectEnd, $nObjectCount, "N", length($asSegments[$i]), "fragment", "yes", ""), "\n";
         }
         else{                       # we have a contig, print an agp contig line and add the sequence to the fsa file
            $nContigCount++;
            my $sContigID = join('_', $sID, $nContigCount);
            print $fhAGP join("\t", $sID, $nObjectStart, $nObjectEnd, $nObjectCount, "W", $sContigID, 1, length($asSegments[$i]), "+"), "\n";
            print $fhFSA ">$sContigID\n$asSegments[$i]\n";
         }
      }
   }
}
