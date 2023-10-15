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
           "version:n"  => \$nGffVersion,
           "prefix:s"   => \$sPrefix);
my $sBedFile = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sBedFile and $sPrefix);
if ($sHelp) {
    my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
    die <<HELP

    $sScriptName [-p -v] <bedfile>

    Takes a bed file and converts it into a gff file of introns and exons
    
    Arguments:
    -v -version
      The GFF version. Default: $nGffVersion
    -p -prefix <sting>
      Output file prefix. Default: $sPrefix

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
      if(scalar(@anBsizes) != $nBcount) { die "Error: mismatch in number of blocks and block sizes\n"; }
      
      # Step through all exon blocks
      for (my $i=0 ; $i<@anBsizes ; $i++){
         
         # Prepare exon coordinates
         my $sExonID = $sStrand eq "+" ? join(":", $i+1, $nBcount) : join(":", $nBcount-$i, $nBcount);
         my $sEname  = $nGffVersion == 1 ? "$sName.ex${sExonID}" : "gene_id \"$sName\"; transcript_id \"$sName\"; exon_id \"$sExonID\";";
         my $nEstart = $nStart + $anBstarts[$i] + 1;
         my $nEend   = $nStart + $anBstarts[$i] + $anBsizes[$i];
         
         # Prepare intron coordinates and print intron line
         if ( $i>0) { 
            my $sIntronID = $sStrand eq "+" ? join(":", $i, $nBcount-1) : join(":", $nBcount-$i, $nBcount-1);
            my $sIname    = $nGffVersion == 1 ? "$sName.in${sIntronID}" : "gene_id \"$sName\"; transcript_id \"$sName\"; intron_id \"$sIntronID\";";
            my $nIstart   = $nStart + $anBstarts[$i-1] + $anBsizes[$i-1] + 1;
            my $nIend     = $nEstart-1;
            print join("\t", $sChr, 'bed2gff', 'intron', $nIstart, $nIend, $nScore, $sStrand, '.', $sIname), "\n";
         }
         
         # Print exon line
         print join("\t", $sChr, 'bed2gff', 'exon', $nEstart, $nEend, $nScore, $sStrand, '.', $sEname), "\n";

      }
   }
   else{ # Partial spec
      die "Error: need bed12 format\n";
   }
}
