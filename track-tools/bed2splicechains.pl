#!/usr/bin/env perl

# Takes a bed file and converts it to a seq of splice chains

use strict;
use warnings;
use Getopt::Long;
use ValidateFiletypes qw(check_bed);

# ARGUMENTS
my $sHelp        = 0;
my $flReplace    = 0;
GetOptions("help!"      => \$sHelp,
           "replace!"   => \$flReplace);
my $sBedFile = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sBedFile);
if ($sHelp) {
    my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
    die <<HELP

    $sScriptName [-r] <bedfile>

    Takes a bed file and converts it into a set of splicechains
    
    Arguments:
    -r -replace
      Replace the bed line name by a unique number

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
my %hSpliceChains;
my ($nMonoCount, $nMultiCount) = (0,0);
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
   $sName = $nCount++ if ($flReplace);
   
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
        
      # Skip mono-exons
      if (scalar @anBsizes > 1){
         my @anSpliceSites;
         for (my $i=0 ; $i<@anBsizes ; $i++){
            my $nEstart = $nStart + $anBstarts[$i];
            my $nEend   = $nStart + $anBstarts[$i] + $anBsizes[$i];
            
            # For the first exon we only keep the exon end
            if ($i == 0){
               push @anSpliceSites, $nEend;
            }
            # For the last exon we only keep the exon start
            elsif ($i == $#anBsizes){
               push @anSpliceSites, $nEstart;
            }
            # For internal exons we keep both ends
            else{
               push @anSpliceSites, $nEstart;
               push @anSpliceSites, $nEend;
            }
         }
         
         # Make the splice chain string and add it to the hash
         my $sSpliceChain = join(":", $sChr, (sort {$a <=> $b} @anSpliceSites), $sStrand);
         $hSpliceChains{$sSpliceChain}{$sName}++;
         $nMultiCount++;
      }
      else{
         $nMonoCount++;
      }
   }
   else{ # Partial spec
      die "Error: need bed12 format\n";
   }
}

# Print the collected splice chains
my $nSpliceChainCount = 0;
foreach my $sSpliceChain (keys %hSpliceChains){
   my $sNames = join(",", sort keys %{$hSpliceChains{$sSpliceChain}});
   print join("\t", $sSpliceChain, $sNames), "\n";
   $nSpliceChainCount++;
}

# Give some final stats
warn ("Found $nSpliceChainCount unique splice chains for $nMultiCount multi-exonic transcripts. Skipped $nMonoCount mono-exonic transcripts.\n") ;
