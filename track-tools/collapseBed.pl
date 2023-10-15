#!/usr/bin/env perl

# Takes a bed file with multiple overlapping transcripts and prepares
# a non-redundant version by a variety of methods, such that each exon annotation
# only occurs in one transcript.

use strict;
use warnings;
use Getopt::Long;
use ValidateFiletypes qw(check_bed);

# GLOBALS
$ENV{SORT}  ||= 'sort';

# ARGUMENTS
my $sHelp        = 0;
my $sMethod      = 'merge';
GetOptions("help!"      => \$sHelp,
           "method:s"   => \$sMethod);

# PRINT HELP
$sHelp = 1 unless($sMethod and @ARGV);
if ($sHelp) {
    my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
    die <<HELP

    $sScriptName [-m] <bedfile1> .. <bedfileN>

    Takes one or more bed files with multiple overlapping transcripts and prepares a 
    non-redundant version such that each exon annotation only occurs in one transcript.
    Note that the thick start and thick end information is lost in the conversion.

    Required:
    -m -method <string>
      The method used for selecting non-redundant transcripts, one of:
       merge    :  Merge exons of overlapping transcripts into a new nr transcript
       largest  :  Pick the largest transcripts based on number of exons & size
      Default: $sMethod

    -help
      This help message
      
HELP
}


#########
# START #
#########

# Check arguments
my $rfGetTranscripts;
$sMethod = lc($sMethod);
if ($sMethod eq 'merge'){
   $rfGetTranscripts = \&get_merged_transcript;
}
elsif ($sMethod eq 'largest'){
   $rfGetTranscripts = \&get_largest_transcript;
}
else{
   die "Error: unknown selection method '$sMethod'\n";
}

# Check the bed file format
foreach my $sBedFile (@ARGV){
   my @asBedErrors = check_bed($sBedFile, 12);
   if (@asBedErrors){
      unshift @asBedErrors, "The following errors were found in bed file '$sBedFile':";
      die join("\n", @asBedErrors), "\n";
   }
}

# Sort the bed files and process the results
my %hLastPos;
my $sInputFiles = join(' ', @ARGV);
open INPUT, "$ENV{SORT} -S 2G -t '\t' -k1,1 -k2,2n -k3,3n $sInputFiles |" or die "Can't sort input: $!\n";
while (<INPUT>){
   next if /^\s*$/;
   next if /^\s*#/;
   s/[\n\r]//g;
   my (@asBedLine) = (split /\t/);
   my ($sChr, $nStart, $nEnd, $sStrand) = @asBedLine[0,1,2,5];
   if(exists($hLastPos{$sStrand})){
      if ( ($sChr ne $hLastPos{$sStrand}{'chr'}) or ($nStart > $hLastPos{$sStrand}{'end'}) ){
         print join("\n", $rfGetTranscripts->($hLastPos{$sStrand}{'transcripts'})), "\n";
         
         $hLastPos{$sStrand}{'chr'}         = $sChr;
         $hLastPos{$sStrand}{'start'}       = $nStart;
         $hLastPos{$sStrand}{'end'}         = $nEnd;
         $hLastPos{$sStrand}{'transcripts'} = [ [@asBedLine] ];  # Initialize first element of array of arrays
      }   
      else{
         $hLastPos{$sStrand}{'end'}  = $nEnd if ($nEnd > $hLastPos{$sStrand}{'end'});
         push @{ $hLastPos{$sStrand}{'transcripts'} }, [@asBedLine];
      }
   }
   else{
      $hLastPos{$sStrand}{'chr'}         = $sChr;
      $hLastPos{$sStrand}{'start'}       = $nStart;
      $hLastPos{$sStrand}{'end'}         = $nEnd;
      $hLastPos{$sStrand}{'transcripts'} = [ [@asBedLine] ];  # Initialize first element of array of arrays
   }
}
# Still need to print the last annotations here!
foreach my $sStrand (keys(%hLastPos)){
   print join("\n", $rfGetTranscripts->($hLastPos{$sStrand}{'transcripts'})), "\n";
}

###############
# SUBROUTINES #
###############

# get_largest_transcript(@aaBedLineArray)
#
# Sorts an array of arrays of bed file lines by the number of blocks
# and returns the transcript with the largest number of blocks
sub get_largest_transcript{
   my $rTranscripts = shift @_;
   my @aaSorted = sort { $b->[9] <=> $a->[9] || $b->[2] <=> $a->[2] } @$rTranscripts;
   return( join("\t", @{$aaSorted[0]}) );
}


# get_merged_transcripts
#
# Merge transcripts with overlapping exons on the same strand together
sub get_merged_transcript{
   my $rTranscripts = shift @_;
   my %hCommonFields;

   # Extract exons with transcript names in array of arrays
   my @aaExons;
   foreach my $rBedLine (@$rTranscripts){
      my ($sChr, $nStart, $nEnd, $sName, $nScore, $sStrand, $nTstart, $nTend, $sRgb, $nBcount, $sBsizes, $sBstarts) = @$rBedLine;
      $hCommonFields{'chr'}    ||= $sChr;
      $hCommonFields{'strand'} ||= $sStrand;
      if($sChr ne $hCommonFields{'chr'})    { die "Error: Features map to different chromosomes in get_merged_transcript\n"; }
      if($sStrand ne $hCommonFields{'strand'}) { die "Error: Features map to different strands ($sStrand, $hCommonFields{'strand'}) in get_merged_transcript\n"; }
      $sBsizes  =~ s/,$//; # Remove trailing comma
      $sBstarts =~ s/,$//; # Remove trailing comma
      my @anBsizes  = split /\s*,\s*/, $sBsizes;
      my @anBstarts = split /\s*,\s*/, $sBstarts;
      if(scalar(@anBsizes) != scalar(@anBstarts)) { die "Error: mismatch in number of block sizes and starts\n"; }
      for (my $i=0 ; $i<@anBsizes ; $i++){
         my $nEstart = $nStart + $anBstarts[$i];
         my $nEend   = $nStart + $anBstarts[$i] + $anBsizes[$i];
         push @aaExons, [$nEstart, $nEend, $sName];
      }
   }
   
   # Sort exons by start position and merge overlapping exons together
   # Keep a hash, linking the exon number to associated transcript IDs
   my @aaMergedExons;
   my %hExonIDtoTranscriptIDs;
   @aaExons = sort{ $a->[0] <=> $b->[0] } @aaExons;
   
   my $rFirstExon = shift @aaExons;
   my $nLastStart = $rFirstExon->[0];
   my $nLastEnd   = $rFirstExon->[1];
   my @asTranscriptIDs  = ($rFirstExon->[2]);
   foreach my $rExon (@aaExons){
      my ($nEstart, $nEend, $sTranscriptID) = @$rExon;
      if ($nEstart >$nLastEnd){
         push @aaMergedExons, [$nLastStart, $nLastEnd];
         $hExonIDtoTranscriptIDs{$#aaMergedExons} = [@asTranscriptIDs];
         $nLastStart = $nEstart;
         $nLastEnd   = $nEend;
         @asTranscriptIDs  = ($sTranscriptID);
      }
      else{
         $nLastEnd   = $nEend if ($nEend > $nLastEnd);
         push @asTranscriptIDs, $sTranscriptID;
      }
   }
   push @aaMergedExons, [$nLastStart, $nLastEnd];
   $hExonIDtoTranscriptIDs{$#aaMergedExons} = [@asTranscriptIDs];
   
   # Group transcripts with overlapping exons together
   my %hTranscriptIDtoGroup;
   my $nGroupCounter = 0;
   foreach my $rTranscriptIDs (values(%hExonIDtoTranscriptIDs)){
      # First check if any of the transcripts have already been assigned to a group before
      my %hExistingGroups;
      foreach my $sTranscriptID (@$rTranscriptIDs){
         $hExistingGroups{$hTranscriptIDtoGroup{$sTranscriptID}}++ if (exists($hTranscriptIDtoGroup{$sTranscriptID}));
      }
      
      # Assign IDs to groups based on exon overlap
      if (scalar(keys(%hExistingGroups)) == 0){     # None of the current IDs has been associated to a group yet, assign new
         foreach my $sTranscriptID (@$rTranscriptIDs){
            $hTranscriptIDtoGroup{$sTranscriptID} = $nGroupCounter;
         }
         $nGroupCounter++;
      }
      elsif (scalar(keys(%hExistingGroups)) == 1){  # One or more IDs are associated to a single group, assign rest to same group
         my $nGroupID = (keys(%hExistingGroups))[0];
         foreach my $sTranscriptID (@$rTranscriptIDs){
            $hTranscriptIDtoGroup{$sTranscriptID} = $nGroupID;
         }
      }
      else{                                         # Multiple groups were found, merge everything into a new group
         while( my ($sTranscriptID, $nGroupID) = each(%hTranscriptIDtoGroup) ){ # Remap existing groups
            $hTranscriptIDtoGroup{$sTranscriptID} = $nGroupCounter if(exists($hExistingGroups{$nGroupID}));
         }
         foreach my $sTranscriptID (@$rTranscriptIDs){ # Map groups for new transcript IDs
            $hTranscriptIDtoGroup{$sTranscriptID} = $nGroupCounter;
         }
         $nGroupCounter++;
      }
   }

   # Assign exons to transcript clusters based on the first transcript ID associated with the exon
   my %hGroupToExons;
   foreach my $nExonID (keys(%hExonIDtoTranscriptIDs)){
      my $nGroupID = $hTranscriptIDtoGroup{$hExonIDtoTranscriptIDs{$nExonID}->[0]};
      push @{ $hGroupToExons{$nGroupID} }, $nExonID;
   }

   # Get a list of transcript IDs belonging to each group
   my %hGroupToTranscriptIDs;
   while( my ($sTranscriptID, $nGroupID) = each(%hTranscriptIDtoGroup) ){
      push @{ $hGroupToTranscriptIDs{$nGroupID} }, $sTranscriptID;
   }

   # Output one line per cluster
   my @asResult;
   foreach my $nGroupID (keys(%hGroupToExons)){
      my @anExonIDs = @{$hGroupToExons{$nGroupID}};
      
      my @aaExonSet = sort { $a->[0] <=> $b->[0] } @aaMergedExons[@anExonIDs];
      my $nBedStart = $aaExonSet[0][0];
      my $nBedEnd   = $aaExonSet[$#aaExonSet][1];
      my $nBlocks   = scalar(@aaExonSet);
      my $sBlockStarts;
      my $sBlockSizes;
      foreach my $rExon (@aaExonSet){
         my ($nEstart, $nEend) = @$rExon;
         $sBlockStarts .= $nEstart-$nBedStart . ',';
         $sBlockSizes  .= $nEend-$nEstart . ',';
      }
      my $sBedID       = join("|", @{$hGroupToTranscriptIDs{$nGroupID}});
      push @asResult, join("\t", $hCommonFields{'chr'}, $nBedStart, $nBedEnd, $sBedID, 0, $hCommonFields{'strand'},
                                 $nBedStart, $nBedEnd, 0, $nBlocks, $sBlockSizes, $sBlockStarts);
   }
   
   return @asResult;
}
