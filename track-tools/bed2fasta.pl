#!/usr/bin/env perl

# Extract DNA sequences from a genome based on supplied bed file
# Note that this script loads the full genome into memory, so it 
# will require at least 3.5 Gb RAM for mammalian genomes.

# GLOBALS
my $nMaxWarnings = 15;  # Maximum number of warnings reported

# MODULES
use strict;
use warnings;
use Getopt::Long;
use ValidateFiletypes qw(check_bed);

# ARGUMENTS
my $sHelp         = 0;
my $sBedFile      = '';
my $sGenomeFile   = '';
my $flFullRegion  = 0;
my $nBasesPerLine = 50;
GetOptions("help!"       => \$sHelp,
           "bedfile:s"   => \$sBedFile,
           "genome:s"    => \$sGenomeFile,
           "fullregion!" => \$flFullRegion,
           "length:i"    => \$nBasesPerLine);

# PRINT HELP
$sHelp = 1 unless($sBedFile and $sGenomeFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

    $sScriptName [-f] -j <bed-file> -g <genome-file>

   Extract bed feature sequences from a genome. By default it will 
   extract exonic sequences only. Add the '-f' parameter to get the 
   sequence for the complete region.

    Required:
    -b <string>
      bed file (full specification)
    -g <string>
      Genome file in fasta format. The fasta headers must contain
      chromosome names that match those used in the bed file.
    -l <integer>
      The number of bases on each line. Default: $nBasesPerLine
    -f
      Output the sequence of the complete bed region (from start to end)
      instead of just the exon sequences.

    -help
      This help message
      
HELP
}


#########
# START #
#########

# Check arguments
die "Error: the number of bases per line must be a positive integer\n" unless ($nBasesPerLine =~ /^[1-9]\d*$/);

# Check whether the bed file is properly formatted
my $nRequiredFields = $flFullRegion ? 4 : 12;
my @asBedErrors = check_bed($sBedFile, $nRequiredFields);
if (@asBedErrors){
   unshift @asBedErrors, "The following errors were found in the input bed file:";
   die join("\n", @asBedErrors), "\n";
}

# Read complete genome sequence into memory
my %hGenomeSeq;
&read_genome($sGenomeFile, \%hGenomeSeq);
print STDERR "Finished reading genome sequences\n";

# Process the bed file and extract the junction sequences
my @asWarnings;
my $nSkippedCount;
my %hUniqueIDcheck;
open BED, $sBedFile or die "Error: can't open '$sBedFile': $!\n";
while (<BED>){
   next if /^\s*#/;
   next if /^track/;
   s/[\n\r]//g;
   my @asBedLine = split /\t/;
   my ($sChr, $nStart, $nEnd, $sName, $nScore, $sStrand, $nTstart, $nTend, $sRgb, $nBcount, $sBsizes, $sBstarts) = @asBedLine;

   # Check for duplicate IDs, no need to check other stuff as it is already handled by the bed file check
   die "Error: the junction ID '$sName' occurs multiple times in the junction file\n" if (exists $hUniqueIDcheck{$sName} );
   $hUniqueIDcheck{$sName}++;

   # Add the junction sequences to the bed output
   if (exists($hGenomeSeq{$sChr})){
      if ($nEnd <= length($hGenomeSeq{$sChr}) ){
         my $sSeq = "";
         if ($flFullRegion){
            $sSeq = substr($hGenomeSeq{$sChr}, $nStart, $nEnd - $nStart);
         }
         else{
            $sBsizes  =~ s/,$//; # Remove trailing comma
            $sBstarts =~ s/,$//; # Remove trailing comma
            my @anBsizes  = split /\s*,\s*/, $sBsizes;
            my @anBstarts = split /\s*,\s*/, $sBstarts;
            die "Error: mismatch between number of block starts and block sizes on line $.\n" unless (scalar(@anBsizes) == scalar(@anBstarts));
            for (my $i=0 ; $i<@anBsizes ; $i++){
               $sSeq .= substr($hGenomeSeq{$sChr}, ($nStart + $anBstarts[$i]), $anBsizes[$i]);
            }
         }
         if ($sStrand){
            if    ($sStrand eq "+"){}
            elsif ($sStrand eq "-"){$sSeq = reverse_complement($sSeq)}
            else{ die "Error: unknown strand type on line $.\n" }
         }
         
         # Break the sequence in shorter chunks
         $sSeq =~ s/.{$nBasesPerLine}/$&\n/sg;
         $sSeq =~ s/\n$//; # Make sure to remove any trailing newlines
         print ">$sName\n$sSeq\n";
      }
      else{
         push @asWarnings, "Warning: feature $sName was skipped because it maps outside the available chromosome sequence";
      }
   }
   else{
      push @asWarnings, "Warning: feature $sName was skipped because there is no sequence data for chromosome '$sChr'";
   }
}
close BED;

# Print warnings
my $nWarningCount = 1;
foreach my $sWarning (@asWarnings){
   if ($nWarningCount < $nMaxWarnings){
      warn "$sWarning\n";
   }
   else{
      warn "-- Further warnings suppressed --\n";
      last;
   }
   $nWarningCount++;
}


###############
# SUBROUTINES #
###############

# read_genome(fasta-filename)
#
# Reads the fasta formatted genome and produces a hash with 
# an entry for each chromosome with the full chromosome sequence
sub read_genome{
   my ($sInput, $rhGenomeSeq) = @_;
   
   my $sChr = '';
   open IN, $sInput or die "Error: can't open $sInput: $!\n";
   while(<IN>){
      next if /^\s*$/;
      next if /^\s*#/;
      s/[\n\r]+$//;
      if (/^\s*>/){
         s/^\s+//; s/ .*$//; s/>//g;
         $sChr = $_;
      }
      else{
         $rhGenomeSeq->{$sChr} .= $_;
      }
   }
   close IN;
}


# rev_comp
#
# Returns the reverse-complement of the supplied sequence
sub reverse_complement{
   my $seq = shift(@_);
   my $rev = reverse $seq;
   $rev =~ tr/ACGTacgt/TGCAtgca/;
   return $rev;
}
