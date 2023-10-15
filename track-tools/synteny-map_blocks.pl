#!/usr/bin/env perl


# Script to map tiling array probe measurements from a syntenic genome
# onto another source genome.


# MODULES
use strict;
use Getopt::Long;

# GET ARGUMENTS
my $sHelp         = 0;
my $sInputAxt     = '';
GetOptions("help!"        => \$sHelp);


# PRINT HELP
$sHelp = 1 unless(@ARGV);
if ($sHelp) {
    die <<HELP

    synteny-map_blocks.pl <axt input file1> .. <axt input fileN>

    Prepares a gff file that can be loaded into IGB to show the synteny
    blocks for another genome. It takes a list of axt-formatted files as
    input. These files can be downloaded from the axtNet dir on the UCSC
    ftp server.

    arguments:
    -h Help

HELP
}

###########
## START ##
###########

my %hChrScores;
my $nChrScore = 40;

print "track name=\"synteny\" useScore=0\n";
foreach my $sFile (@ARGV){
   open INPUT, $sFile or die "Can't open .axt file: $!\n";
   while (<INPUT>){
      next if (/^\s*$/);
      next if (/^#/);
      chomp;
      if (/^\d+/){
         my ($nBlockID, $sRefChr, $nRefStart, $nRefEnd, $sSyntChr, $nSyntStart, $nSyntEnd, $sStrand, $nScore) = split /\s/;
         
         # Give each chromosome a different score, and thus a different color in IGB
         # Score increments by 40, which should give space for 25 chromosomes (1000 is max score for gff file)
         unless (exists($hChrScores{$sSyntChr})){
            $hChrScores{$sSyntChr} = $nChrScore;
            $nChrScore += 40;
         }
         print join("\t", $sRefChr, 'axnet', 'synteny', $nRefStart, $nRefEnd, $hChrScores{$sSyntChr}, $sStrand, '.', "$sSyntChr.$nBlockID"), "\n";
      }
   }
}
