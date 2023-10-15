#!/usr/bin/env perl

# Takes a bed file and converts it to a gff file, keeping exon information

use strict;
use warnings;
use Getopt::Long;
use ValidateFiletypes qw(check_bed);

# ARGUMENTS
my $sHelp        = 0;
my $sReplace     = 0;
my $nFlankSize   = 2500;
my $nMinExonSize = 100;
my $nGffVersion  = 3;
my $sPrefix      = "output";
GetOptions("help!"      => \$sHelp,
           "replace!"   => \$sReplace,
           "flank:n"    => \$nFlankSize,
           "version:n"  => \$nGffVersion,
           "prefix:s"   => \$sPrefix,
           "minexon:n"  => \$nMinExonSize);
my $sBedFile = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sBedFile and $sPrefix);
if ($sHelp) {
    my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
    die <<HELP

    $sScriptName [-r -e -v] <bedfile>

    Takes a bed file and converts it into a gff file
    -f -flank <integer>
      Intron flank size. Default: $nFlankSize;
    -r -replace
      Replace the bed line name by a unique number
    -v -version
      The GFF version. Default: $nGffVersion
    -p -prefix <sting>
      Output file prefix. Default: $sPrefix
    -m -minexon <integer>
      Minimum exon size. Default: $nMinExonSize

    -help
      This help message
      
HELP
}


#########
# START #
#########

# Check the bed file
my @asBedErrors = check_bed($sBedFile);
if (@asBedErrors){
   unshift @asBedErrors, "The following errors were found in the input bed file:";
   die join("\n", @asBedErrors), "\n";
}

open EXON, ">${sPrefix}_exon.gff" or die "Error: can't write to ${sPrefix}_exon.gff: $!\n";
open INTRON, ">${sPrefix}_intron.gff" or die "Error: can't write to ${sPrefix}_intron.gff: $!\n";

# Open the bed file
my %hNames;
my $nCount = 1;
open INPUT, $sBedFile or die "Can't open $sBedFile: $!\n";
while (<INPUT>){
   next if /^\s*$/;
   next if /^\s*#/;
   next if /^\s*[Tt]rack/;
   s/[\n\r]//g;
   my @asBedLine = split /\t/;
   my ($sChr, $nStart, $nEnd, $sName, $nScore, $sStrand, $nTstart, $nTend, $sRgb, $nBcount, $sBsizes, $sBstarts) = @asBedLine;
   
   # Warn if an identifier occurs multiple times in the bed file because it can cause problems for gff
   warn "Warning: identifier '$sName' occurs multiple times in the bed file\n" if (exists($hNames{$sName}));
   $hNames{$sName}++;
   $sName = $nCount++ if ($sReplace);
   
   # Set defaults for some of the fields
   $nScore  = 0 unless($nScore);
   $sStrand ||= '.';
   $nTstart ||= $nStart;
   $nTend   ||= $nEnd;
   
   # Distinguish between a full and a partial bed file specification
   if ($sBsizes and $sBstarts){ # Full spec
      $sBsizes  =~ s/,$//; # Remove trailing comma
      $sBstarts =~ s/,$//; # Remove trailing comma
      my @anBsizes  = split /\s*,\s*/, $sBsizes;
      my @anBstarts = split /\s*,\s*/, $sBstarts;
      if(scalar(@anBsizes) != scalar(@anBstarts)) { die "Error: mismatch in number of block sizes and starts\n"; }
     
      for (my $i=0 ; $i<@anBsizes ; $i++){
         if ($sStrand eq "+"){
            if ($i < $#anBsizes) { # Skip last exon
               my $nExonID = $i + 1;
               my $nEstart = $nStart + $anBstarts[$i];
               my $nEend   = $nStart + $anBstarts[$i] + $anBsizes[$i];
               my $nNextStart = $nStart + $anBstarts[$i+1];
               my $nIntronLen = $nNextStart - $nEend;
               if ( ($nIntronLen >= $nFlankSize) and ( ($nEend-$nEstart) >= $nMinExonSize ) ){
                  my $sOut = "$sName.${nExonID}";
                  $sOut    = "gene_id \"$sOut\"; transcript_id \"$sOut\";" if ($nGffVersion>1);
                  print EXON join("\t", $sChr, 'bed2gff', 'exon', $nEstart+1, $nEend, $nScore, $sStrand, '.', $sOut), "\n";
                  print INTRON join("\t", $sChr, 'bed2gff', 'exon', $nEend, $nEend+$nFlankSize, $nScore, $sStrand, '.', $sOut), "\n";
               }
            }
         }
         else{
            if ($i > 0) { # Skip last exon
               my $nExonID = $nBcount - $i;
               my $nEstart = $nStart + $anBstarts[$i];
               my $nEend   = $nStart + $anBstarts[$i] + $anBsizes[$i];
               my $nPreviousEnd   = $nStart + $anBstarts[$i-1] + $anBsizes[$i-1];
               my $nIntronLen     = $nEstart - $nPreviousEnd;
               if ( ($nIntronLen >= $nFlankSize) and ( ($nEend-$nEstart) >= $nMinExonSize ) ){
                  my $sOut   = "$sName.${nExonID}";
                  $sOut   = "gene_id \"$sOut\"; transcript_id \"$sOut\";" if ($nGffVersion>1);
                  print EXON join("\t", $sChr, 'bed2gff', 'exon', $nEstart+1, $nEend, $nScore, $sStrand, '.', $sOut), "\n";
                  print INTRON join("\t", $sChr, 'bed2gff', 'exon', $nEstart-$nFlankSize, $nEstart, $nScore, $sStrand, '.', $sOut), "\n";
               }
            }
         }
      }
   }
   else{ # Partial spec
      die "Error: need bed12 format\n";
   }
}

close EXON;
close INTRON;
