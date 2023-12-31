#!/usr/bin/env perl

# Parse the blat output files generated by the script map-array-probes_run-blat.pl
# Please note that the array probe identifiers in the blat output file need to be 
# specifically formatted. Please use the map-array-probes_bpmap2tab.pl script to ensure 
# proper formatting of the oligo identifier in the blat output file.

# MODULES
use strict;
use Getopt::Long;

$ENV{'BPMAP2TPMAP_BINARY'} ||= 'bpmap2tpmap';
$ENV{'TPMAP2BPMAP_BINARY'} ||= 'tpmap2bpmap';

# GET ARGUMENTS
my $sHelp           = 0;
my $sStat           = 0;
my $sRestrict       = 0;
my $sPreserveStrand = 'b';
my $sInputBlatDir   = '';
my $sInputBpmap     = '';
my $sOutput         = '';
my $sDbVersion      = '';
my $nThreshold      = 0;
GetOptions("help!"        => \$sHelp,
           "restrict!"    => \$sRestrict,
           "strand:s"     => \$sPreserveStrand,
           "input:s"      => \$sInputBlatDir,
           "bpmap:s"      => \$sInputBpmap,
	   "output:s"     => \$sOutput,
	   "version:s"    => \$sDbVersion,
	   "threshold:s"  => \$nThreshold);


# PRINT HELP
$sHelp = 1 unless($sInputBlatDir and $sOutput);
if ($sHelp) {
    die <<HELP

    map-array-probes_blat2bpmap.pl -b <bpmap input> -i <blat input dir> -o <output name> 

    arguments:
    -b Bpmap input file
        Name of the bpmap file to remap
    -i Blat input directory
        Directory with blat output files (without header). These blat output files need to
        be generated by first running map-array-probes_bpmap2tab.pl (generates a tab-delimited
        file with probe sequences from a bpmap file) and then map-array-probes_run-blat.pl
        (to BLAT the probe sequences against the target genome).
    -o Output file prefix
        Name of the remapped bpmap output file
    -v Genome database version
        The NCBI version number of the genome that the probes are being mapped to, e.g.  NCBIv33
    -r Restrict output
        Restricts output to only include probes that map to the chromosomes the array was 
        originally designed for. All other probes are mapped normally, but excluded as tiling
        array probes at the sequence_group and probe_type level.
    -s Strand information
        Indicate what to do with the strand information for mapped probes:
         b   Keep both strands
         u   Keep both strands but set strand information to 't' for all probes
         f   Keep only reverse strand matches
         r   Keep only forward strand matches
    -t Score threshold
        Threshold for alignment scores. Currently there are 3 levels: 
           1    - Unique perfect match; 
           0.5  - Match selected by neighborhood mapping from multiple perfect matches
           0.2  - Match selected by best guess estimation
    -h Help

HELP
}

###########
## START ##
###########

# Add versioning to output file name
$sOutput = join("_", $sOutput, $sDbVersion);

# Read the blat hits
print "Reading blat hits\n";
my ($rLUT, $sSeqGroup) = &read_blat_hits($sInputBlatDir);

# Spit out some statistics if asked for
if ($sStat){
   print "Processing statistics\n";
   &print_hit_statistics($rLUT);
}

# Read bpmap file content and immediately process unique perfect-match hits
print "Processing bpmap input file\n";
my ($rMappings, $rDuplicatesIdx, $nProbeTot, $nProbeUnk, $nProbePerf, $rOriginalMappedChromosomes) = &read_bpmap($sInputBpmap, $rLUT, $sSeqGroup, $sOutput);
my $nProbeMatch = $nProbeTot - $nProbeUnk;
print "$nProbeTot probes loaded for mapping\n";
print "$nProbeUnk probes without a perfect match in the genome -> moved to 'Unknown' feature group\n";
print "$nProbeMatch probes with at least one single perfect match in the genome -> proceeding with multi-stage mapping\n";
print "Pass 1, mapping perfect-match probes with single genomic match\n";
print join(" ", $nProbePerf, "probes matched in pass 1" ), "\n";

# Process duplicates entries, forward and reverse screen
# For the ones that still cannot be mapped, do best-guess mapping
print join(" ", "Pass 2, forward neighbor-based mapping of probes with multiple genomic matches." ), "\n";
my ($rDuplicatesIdx, $nMatched) = &map_duplicates($rMappings, $rDuplicatesIdx, $rLUT, 'forward');
print join(" ", $nMatched, "probes matched in pass 2" ), "\n";
   
print join(" ", "Pass 3, reverse neighbor-based mapping of probes with multiple genomic matches." ), "\n";
my ($rDuplicatesIdx, $nMatched) = &map_duplicates($rMappings, $rDuplicatesIdx, $rLUT, 'reverse');
print join(" ", $nMatched, "probes matched in pass 3" ), "\n";
  
print join(" ", "Pass 4, best-guess mapping of remaining probes with multiple genomic matches." ), "\n";
my $nMatched = &best_guess_duplicates($rMappings, $rDuplicatesIdx, $rLUT);
print join(" ", $nMatched, "probes matched in pass 4" ), "\n";



# Write remapped bpmap entries to tpmap file
print "Adding remapped entries to temporary '$sOutput.tpmap' file\n";
&tpmap_write_remapped_entries($rMappings, $sSeqGroup, $sDbVersion, $sOutput, $rOriginalMappedChromosomes, $sRestrict, $sPreserveStrand, $nThreshold);

# Convert tpmap to bpmap file
#print "Writing bpmap file: '$sOutput.bpmap'\n";
#`$ENV{'TPMAP2BPMAP_BINARY'} -tpmap $sOutput.tpmap -bpmap $sOutput.bpmap`;


#################
## SUBROUTINES ##
#################

# read_blat_hits($sInputBlatDir)
#
# Processes all blat output files in the specified dir and returns
# a hash of arrays with pm hits for each array probe
sub read_blat_hits {
   my ($sInputBlatDir) = @_;
   
   my %hSeqGroup;
   my %hLUT;
   opendir DIR, $sInputBlatDir or die "Can't open blat output dir: $!\n";
   while (defined(my $sFile = readdir(DIR))) {
      next unless ($sFile =~ /\.psl/);
      open INPUT, "$sInputBlatDir/$sFile" or die "Can't open BLAT output file '$sFile': $!\n";
      while (<INPUT>){
         next if (/^\s*$/);
         s/[\n\r]//g;
         s/ //g;
         my (@asTmp) = split /\t/;
         my ($nMatch, $nGaps, $sStrand, $sQname, $nQsize, $sTname, $nTstart) = @asTmp[0,6,8,9,10,13,15];
         my ($sSequenceGroup, @asIdentifier) = split /\|/, $sQname;
         my $sOligoID       = join('|', @asIdentifier);
         $hSeqGroup{$sSequenceGroup}++;
         if ( ($nMatch == $nQsize) and ($nGaps==0) ){
            push @{$hLUT{$sOligoID}}, [$sTname, $nTstart, $sStrand];
         }
      }
   }
   close INPUT;
      
   my @asSeqGroups = keys(%hSeqGroup);
   die "Multiple sequence groups found in BLAT output, cannot continue remapping\n"  unless(scalar(@asSeqGroups)==1);
   return (\%hLUT, $asSeqGroups[0]);
}



# print_hit_statistics($rLUT)
#
# Prints a breakdown of which kinds of combinations of numbers of genomic
# matches and mismatches exist for all array probes
sub print_hit_statistics {
   my $rLUT = shift @_;
   
   my %hMatches;
   foreach my $sOligoID ( keys(%$rLUT) ){
      $hMatches{scalar(@{$rLUT->{$sOligoID}})}++;
   }
   
   print "# of PM\tCount\n";
   foreach my $nPM ( sort {$a <=> $b} keys(%hMatches)){
         print join ("\t", $nPM, $hMatches{$nPM}), "\n";
   }
}



# read_bpmap($sInputBpmap, $rLUT, $sSeqGroup $sOutput)
#
# Reads the bpmap file and channels the sequence groups that do not need remapping
# directly to the temporary tpmap output file. The sequence groups that do need remapping 
# are read into a single array structure and reliable PM probes are immediately remapped 
# on the fly.
sub read_bpmap{
   my ($sInputBpmap, $rLUT, $sSeqGroup, $sOutput) = @_;
   
   my %hOriginalMappedChromosomes;      # We need to remember which chromosomes were originally mapped
   my ($nProbeTot, $nProbeUnk, $nProbePerf) = (0,0,0);
   my @aaMapping;
   my @anDuplicatesIdx;
   my $flPrint    = 1;
   
   # Open and process the bpmap file
   open OUTPUT, ">$sOutput.tpmap" or die "Can't open temporary output file: $!\n";
   open BPMAP, "$ENV{'BPMAP2TPMAP_BINARY'} -bpmap $sInputBpmap |" or die "Can't extract from $sInputBpmap: $!\n";
   while (<BPMAP>){
      # Check for sequence groups
      if (/^#seq_group_name\t(.*)$/){
         if (lc($1) eq lc($sSeqGroup)){ $flPrint = 0 }
         else                         { $flPrint = 1 }
      }
      
      # Either print to output or read into array
      if ($flPrint){
         print OUTPUT;
      }
      else{
         next if (/^#/);
         next if (/^\s*$/);
         chomp;
         $nProbeTot++;
         
         # Collect previous probe mapping
         my %hProbeMapData;
         my (@asLine) = split /\t/;
         if    (@asLine == 7) { @hProbeMapData{qw(sequence old_strand old_chr old_pos pm_x pm_y score)}           = @asLine }
         elsif (@asLine == 9) { @hProbeMapData{qw(sequence old_strand old_chr old_pos pm_x pm_y mm_x mm_y score)} = @asLine }
         else                 { die "Incorrect number of fields in bpmap file\n"}
         $hOriginalMappedChromosomes{lc($asLine[2])}++;    # We need to remember which chromosomes were originally mapped
         
         # Retrieve mapping data if there is a single unique match (or no match at all)
         my $sID = join('|', $hProbeMapData{pm_x}, $hProbeMapData{pm_y});
         my (@asUniqueMapping) = get_unique_blat_mapping($sID, $rLUT);
         if (@asUniqueMapping){
            @hProbeMapData{qw(new_chr new_pos new_strand score)} = @asUniqueMapping;
            if ( lc($asUniqueMapping[0]) eq 'unknown'){  $nProbeUnk++; }
            else{ $nProbePerf++; }
            push @aaMapping, {%hProbeMapData};
         }
         else{
            push @aaMapping, {%hProbeMapData};
            push @anDuplicatesIdx, $#aaMapping;
         }
         
      }      
   }
   
   close OUTPUT;
   close BPMAP;
   
   return (\@aaMapping, \@anDuplicatesIdx, $nProbeTot, $nProbeUnk, $nProbePerf, \%hOriginalMappedChromosomes);
}



# get_unique_blat_mapping($sID, $rLUT)
# 
# Gets the unique blat mapping for a specified probe identifier as an array
# with elements: chromosome, chromosome_pos, mapping_score
# If more than one perfect-match is found in the genome, an empty array is returned
sub get_unique_blat_mapping{
   my ($sID, $rLUT) = @_;
   my @asResult;
   
   if (exists($rLUT->{$sID})){
      @asResult = (@{$rLUT->{$sID}->[0]}, 1) if ( scalar(@{$rLUT->{$sID}})==1 );
   }
   else{
      @asResult = ('Unknown', 0, '+', 0);
   }
   
   return @asResult;
}



# map_duplicates($rMappings, $rDuplicates, $rLUT, 1)
#
# Attempts to map probes with multiple genomic matches by selecting
# the genomic copy that is closest to neighboring probes in the original mapping
sub map_duplicates {
   my ($rMappings, $rDuplicates, $rLUT, $sDirection) = @_;
   my @anReturn;
   my $nSearchWindow   = 10;   # How many neighboring probes to consider in the search
   my $nDistanceThresh = 200;  # How much extra difference in interprobe distance can there be
   
   # Set search direction
   my @anDuplicates = @$rDuplicates;
   @anDuplicates    = reverse(@anDuplicates) if ($sDirection eq 'reverse');
   
   # Process all duplicate probes
   my $nMatched = 0;
   foreach my $nPid (@anDuplicates){
      # Get original probe ID and position      
      my $sProbeIdentifier = join('|', $rMappings->[$nPid]{pm_x}, $rMappings->[$nPid]{pm_y});
      my $sOldProbeChr     = $rMappings->[$nPid]{old_chr};
      my $nOldProbePos     = $rMappings->[$nPid]{old_pos};
   
      # Define search window
      my $nWindowStart = $nPid - 10;
      my $nWindowStop  = $nPid + 10;
      $nWindowStart = 0 if ($nWindowStart<0);
      $nWindowStop  = scalar(@$rMappings) if ($nWindowStop > scalar(@$rMappings));
      
      # Process all neighboring probes in the window
      my ($flRemapped, $nRemapID, $nRemapDistDiff) = (0,0,0);
      for (my $i=$nWindowStart ; $i<$nWindowStop ; $i++){
         next unless (exists($rMappings->[$i]{new_chr}));               # neighbor was not mapped
         next if (lc($rMappings->[$i]{new_chr}) eq 'unknown');          # neighbor could not be found in genome
         next if (lc($rMappings->[$i]{new_chr}) ne lc($sOldProbeChr));  # neighbor was originally on a different chromosome
         
         # Extract neighboring probe positions
         my ($nNeighOldProbeChr, $nNeighOldProbePos) = ($rMappings->[$i]{old_chr}, $rMappings->[$i]{old_pos});
         my ($nNeighNewProbeChr, $nNeighNewProbePos) = ($rMappings->[$i]{new_chr}, $rMappings->[$i]{new_pos});
         my $nOldProbeDistance = abs($nOldProbePos - $nNeighOldProbePos);

         # Now process all duplicate genomic matches, pick the neighborhood match with the closest match to 
         # the original interprobe distance, not just the first match within the threshold!
         for (my $j=0 ; $j < @{$rLUT->{$sProbeIdentifier}} ; $j++ ){
            my ($sMatchChr, $nMatchPos, $sMatchStrand) = @{$rLUT->{$sProbeIdentifier}[$j]};
            my $nNewProbeDistance = abs($nNeighNewProbePos - $nMatchPos);
            
            if ( ($sMatchChr eq $nNeighNewProbeChr) and ($nNewProbeDistance <= ($nOldProbeDistance+$nDistanceThresh)) ){
               my $nTmpDistDiff = abs($nOldProbeDistance - $nNewProbeDistance);
               if ($flRemapped){
                  if ($nTmpDistDiff<$nRemapDistDiff){
                     $nRemapID       = $j;
                     $nRemapDistDiff = $nTmpDistDiff;
                  }
               }
               else{
                  $nRemapID       = $j;
                  $nRemapDistDiff = $nTmpDistDiff;
                  $flRemapped     = 1;
               }
            }
         }
      }
      
      # Now add the mapping info to the mapping array
      if ($flRemapped){
         my ($sMatchChr, $nMatchPos, $nMatchStrand) = @{$rLUT->{$sProbeIdentifier}[$nRemapID]};
         $rMappings->[$nPid]{new_chr}    = $sMatchChr;
         $rMappings->[$nPid]{new_pos}    = $nMatchPos;
         $rMappings->[$nPid]{new_strand} = $nMatchStrand;
         $rMappings->[$nPid]{score}      = 0.5;
         $nMatched++;
      }
      else{
         push @anReturn, $nPid;
      }
   }
   
   # Return the probes that could not be remapped this turn
   return (\@anReturn, $nMatched);
}



# best_guess_duplicates($rMappings, $rDuplicatesIdx, $rLUT)
#
# Last-attempt mapping of probes that could not be mapped relative to
# neighboring probes in the original mapping. Mapping sequence is as follows:
#  1 - Try to map the probe closest to the position in the original mapping 
#      (i.e. chromosome and position)
#  2 - If this fails, select the first of the duplicate probe mappings
sub best_guess_duplicates {
   my ($rMappings, $rDuplicatesIdx, $rLUT) = @_;
   my $nMatched = 0;

   foreach my $nPid (@$rDuplicatesIdx){
      # Get original probe ID and position
      my $sProbeIdentifier = join('|', $rMappings->[$nPid]{pm_x}, $rMappings->[$nPid]{pm_y});
      my $sOldProbeChr     = $rMappings->[$nPid]{old_chr};
      my $nOldProbePos     = $rMappings->[$nPid]{old_pos};
      
      # Process all the remaining duplicate mappings and try to select the closest match on same chromosome as in original mapping
      my ($sRemapChr, $nRemapPos, $sRemapStrand) = ('', 0, '');
      foreach my $rMatch (@{$rLUT->{$sProbeIdentifier}}){
         my ($sMatchChr, $nMatchPos, $sMatchStrand) = @$rMatch;
         if ( lc($sMatchChr) eq lc($sOldProbeChr) ){
            if ($nRemapPos){
               $nRemapPos = $nMatchPos if (abs($nMatchPos-$nOldProbePos) <= abs($nRemapPos-$nOldProbePos));
            }
            else{
               ($sRemapChr, $nRemapPos, $sRemapStrand) = ($sMatchChr, $nMatchPos, $sMatchStrand);
            }
         }
      }
      
      # Select first mapping if no match could be found on original chromosome
      unless ($sRemapChr){
         ($sRemapChr, $nRemapPos, $sRemapStrand) = @{$rLUT->{$sProbeIdentifier}[0]};
      }
      
      # Add mapping to the output
      $rMappings->[$nPid]{new_chr}    = $sRemapChr;
      $rMappings->[$nPid]{new_pos}    = $nRemapPos;
      $rMappings->[$nPid]{new_strand} = $sRemapStrand;
      $rMappings->[$nPid]{score}      = 0.2;
      $nMatched++;
   }
   
   return ($nMatched);
}



# tpmap_write_remapped_entries($rMappings, $sSeqGroup, $sDbVersion, $sOutput, $rOriginalMappedChromosomes, $sRestrict)
#
# Appends the remapped entries to the tpmap file
sub tpmap_write_remapped_entries {
   ($rMappings, $sSeqGroup, $sDbVersion, $sOutput, $rOriginalMappedChromosomes, $sRestrict, $sPreserveStrand, $nThreshold) = @_;

   # Sort the mapping array, ascending by new chromosome (lowercase, lexographic) and position (numerical)
   my @aaSortedMappings = sort { lc $a->{new_chr} cmp lc $b->{new_chr} || $a->{new_pos} <=> $b->{new_pos} } @$rMappings;
   
   # Append the sorted array to the tpmap file, make a separate sequence group for 'Unknown'
   my $sChr = '';
   open OUTPUT, ">>$sOutput.tpmap" or die "Can't open temporary output file: $!\n";
   foreach my $rMapping (@aaSortedMappings){
      my %hProbeMapData = %$rMapping;
      
      # Set sequence group name and probeset type
      my $sSeqGroupName = $sSeqGroup;
      my $sProbeSetType = 'tiling';
      if ($hProbeMapData{new_chr} eq 'Unknown'){
         $sSeqGroupName = 'Unknown';
         $sProbeSetType = 'arbitrary';
      }
      
      # Optionally disregard probes mapping to chromosomes that were not part of the original array layout
      if ($sRestrict){
         unless(exists($rOriginalMappedChromosomes->{lc($hProbeMapData{new_chr})})) {
            $sSeqGroupName = 'Unknown';
            $sProbeSetType = 'arbitrary';
         }
      }
      
      # Print output
      if (lc($hProbeMapData{new_chr}) ne lc($sChr)){
         print OUTPUT "#seq_group_name\t$sSeqGroupName\n#version\t$sDbVersion\n#probeset_type\t$sProbeSetType\n";
         $sChr = $hProbeMapData{new_chr};
      }
      
      # Recode strand info and check if strandedness actually needs to be preserved
      if    ($hProbeMapData{new_strand} eq '+') { $hProbeMapData{new_strand} = 't'; }
      elsif ($hProbeMapData{new_strand} eq '-') { $hProbeMapData{new_strand} = 'f'; }
      else { die "Unknown strand type `$hProbeMapData{new_strand}` for probe with PM x/y position: $hProbeMapData{pm_x} / $hProbeMapData{pm_y}\n" }
      
      
      # Figure out which lines to print based on strand stuff
      my $nPrintLine = 0;
      if    (lc($sPreserveStrand) eq 'b'){ $nPrintLine = 1; }
      elsif (lc($sPreserveStrand) eq 'u'){ $hProbeMapData{new_strand} = 't'; $nPrintLine = 1; }
      elsif (lc($sPreserveStrand) eq 'f'){ $nPrintLine = 1 if ($hProbeMapData{new_strand} eq 'f'); }
      elsif (lc($sPreserveStrand) eq 'r'){ $nPrintLine = 1 if ($hProbeMapData{new_strand} eq 't'); }
      else{ die "Unknown option: -s '$sPreserveStrand'\n"; }
      
      
      # Print output, depending on whether MM probes are present and whether threshold is met
      if ($nPrintLine and ($hProbeMapData{score} >= $nThreshold) ){
         if (exists($hProbeMapData{mm_x})){
            print OUTPUT join("\t", @hProbeMapData{qw(sequence new_strand new_chr new_pos pm_x pm_y mm_x mm_y score)}), "\n";
         }
         else{
            print OUTPUT join("\t", @hProbeMapData{qw(sequence new_strand new_chr new_pos pm_x pm_y score)}), "\n";
         }
      }
      
   }
   close OUTPUT;
   
}





