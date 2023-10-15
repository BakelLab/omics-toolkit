#!/usr/bin/env perl

# 16.08.2011 14:38:52 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $sReference   = 'query';
GetOptions("help!"       => \$sHelp,
           "reference:s" => \$sReference);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName <psl file> ... <psl file N>
   
   Gives information on the proportion of each query sequence
   that is covered by blat hits.
   
   Options:
    -r --reference <string>
      Set the reference for the coverage analysis ('query' or 'target')
      default: $sReference
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check arguments
$sReference = lc($sReference);
die "Error: unknown reference selection. Choose either 'query' or 'target'\n" unless( $sReference =~ /^(query|target)$/);

# Collect the required info from the psl files
my %hHits;
foreach my $sInput (@ARGV){
   open IN, $sInput or die "Error: can't open '$sInput': $!\n";
   while (<IN>){
      next unless /^\d+/;
      s/[\n\r]+$//;
      my($nMatch,$nMismatch,$nRepMatch,$nNs,$nQgapCount,$nQgapBases,$nTgapCount,$nTgapBases,$sStrand,$sQname,$nQsize,$nQstart,$nQend,
         $sTname,$nTsize,$nTstart,$nTend, $nBlockCount,$sBlockSizes,$sQstarts,$sTstarts) = split /\t/;
      
      my ($sRname,$nRsize,$nRstart,$nRend) = $sReference eq 'query' ? ($sQname,$nQsize,$nQstart,$nQend) : ($sTname,$nTsize,$nTstart,$nTend);
      if ($sRname and $nRsize =~ /^\d+$/ and $nRstart =~ /^\d+$/ and $nRend =~ /^\d+$/){
         push @{$hHits{$sRname}{hits}}, [$nRstart, $nRend];
         if (exists $hHits{$sRname}{size}){
            die "Error: size mismatch for '$sRname' ($nRsize vs. $hHits{$sRname}{size}) on line $.\n" unless ($hHits{$sRname}{size} == $nRsize);
         }
         else{
            $hHits{$sRname}{rsize} = $nRsize;
         }
      }
      else{
         die "Error: line $. in '$sInput' does not conform to the psl format\n";
      }
   }
   close IN;
}

# Now process the hits into segments
my %hSegments;
foreach my $sRname (keys %hHits){
   my @anHits = @{$hHits{$sRname}{hits}};
   @anHits = sort {$a->[0] <=> $b->[0]} @anHits;
   $hSegments{$sRname}{rsize} = $hHits{$sRname}{rsize};
   
   # Merge overlapping hits
   my ($nSegStart, $nSegEnd) = @{shift @anHits};
   push @{$hSegments{$sRname}{segments_missed}},  join('-', 0, $nSegStart-1) if ($nSegStart);
   foreach my $rHit (@anHits){
      my ($nHitStart, $nHitEnd) = @$rHit;
      if ($nHitStart > $nSegEnd){
         push @{$hSegments{$sRname}{segments_matched}}, join('-', $nSegStart, $nSegEnd);
         push @{$hSegments{$sRname}{segments_missed}},  join('-', $nSegEnd+1, $nHitStart-1);
         $hSegments{$sRname}{segsize} += $nSegEnd - $nSegStart;
         ($nSegStart, $nSegEnd) = ($nHitStart, $nHitEnd);
      }
      else{
         $nSegEnd = $nHitEnd > $nSegEnd ? $nHitEnd : $nSegEnd;
      }
   }
   push @{$hSegments{$sRname}{segments_matched}}, join('-', $nSegStart, $nSegEnd);
   if ($nSegEnd < $hSegments{$sRname}{rsize}){
      push @{$hSegments{$sRname}{segments_missed}}, join('-', $nSegEnd+1, $hSegments{$sRname}{rsize});
   }
   $hSegments{$sRname}{segsize} += $nSegEnd - $nSegStart;
}

# And finally output the segment matches
foreach my $sRname (keys %hSegments){
   my $nFraction = sprintf("%0.2f", $hSegments{$sRname}{segsize} / $hSegments{$sRname}{rsize});
   my $sMatched  = join(',', @{$hSegments{$sRname}{segments_matched}});
   my $sMissed   = exists($hSegments{$sRname}{segments_missed}) ? join(',', @{$hSegments{$sRname}{segments_missed}}) : '';
   print join("\t", $sRname, $hSegments{$sRname}{rsize}, $hSegments{$sRname}{segsize}, $nFraction, $sMatched, $sMissed), "\n";
}
