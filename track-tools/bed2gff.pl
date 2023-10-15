#!/usr/bin/env perl

# Takes a bed file and converts it to a gff file, keeping exon information

use strict;
use warnings;
use Getopt::Long;
use ValidateFiletypes qw(check_bed);

# ARGUMENTS
my $sHelp       = 0;
my $sReplace    = 0;
my $flExonOnly  = 0;
my $nGffVersion = 1;
GetOptions("help!"      => \$sHelp,
           "replace!"   => \$sReplace,
           "exon_only!" => \$flExonOnly,
           "version:n"  => \$nGffVersion);
my $sBedFile = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sBedFile);
if ($sHelp) {
    my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
    print <<HELP;

    $sScriptName [-r -e -v] <bedfile>

    Takes a bed file and converts it into a gff file

    -r -replace
      Replace the bed line name by a unique number
    -e -exon_only
      Only print exons, skip CDS entries
    -v -version
      The GFF version
      default: $nGffVersion

    -help
      This help message
      
HELP
exit 0;
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

   # Format identifier depending on the gff version type
   $sName = "gene_id \"$sName\"; transcript_id \"$sName\";" if ($nGffVersion>1);
   
   # Distinguish between a full and a partial bed file specification
   if ($sBsizes and $sBstarts){ # Full spec
      $sBsizes  =~ s/,$//; # Remove trailing comma
      $sBstarts =~ s/,$//; # Remove trailing comma
      my @anBsizes  = split /\s*,\s*/, $sBsizes;
      my @anBstarts = split /\s*,\s*/, $sBstarts;
      if(scalar(@anBsizes) != scalar(@anBstarts)) { die "Error: mismatch in number of block sizes and starts\n"; }
      
      # Process the exon blocks
      my $nCDSsize = 0;
      
      for (my $i=0 ; $i<@anBsizes ; $i++){
         my $nEstart = $nStart + $anBstarts[$i];
         my $nEend   = $nStart + $anBstarts[$i] + $anBsizes[$i]; 
         print join("\t", $sChr, 'bed2gff', 'exon', $nEstart+1, $nEend, $nScore, $sStrand, '.', $sName), "\n";
         
         # Print CDS entries if requested and if the thick start and thick end is different from feature start and end
         unless( (($nTstart==$nStart) and ($nTend==$nEnd)) or $flExonOnly ){
            unless( ($nEend < $nTstart) or ($nEstart > $nTend) ){
               my $nCDSstart = $nTstart > $nEstart ? $nTstart : $nEstart;
               my $nCDSend   = $nTend   < $nEend   ? $nTend   : $nEend;
               my $nFrame    = $nCDSsize % 3;
               $nFrame = $nFrame ? (2-$nFrame)+1 : $nFrame;
               print join("\t", $sChr, 'bed2gff', 'start_codon', $nCDSstart+1, $nCDSstart+3, $nScore, $sStrand, '.', $sName), "\n" if ($nCDSstart == $nTstart);
               print join("\t", $sChr, 'bed2gff', 'CDS', $nCDSstart+1, $nCDSend, $nScore, $sStrand, $nFrame, $sName), "\n";
               $nCDSsize += $nCDSend - $nCDSstart;
            }
         }
      }
   }
   else{ # Partial spec
      print join("\t", $sChr, 'bed2gff', 'exon', $nStart+1, $nEnd, $nScore, $sStrand, '.', $sName), "\n";
      
      # Print CDS entries if requested and if the thick start and thick end is different from feature start and end
      unless( (($nTstart==$nStart) and ($nTend==$nEnd)) or $flExonOnly ){
         print join("\t", $sChr, 'bed2gff', 'CDS', $nTstart+1, $nTend, $nScore, $sStrand, 0, $sName), "\n";
      }
   }
}
