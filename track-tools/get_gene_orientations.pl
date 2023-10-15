#!/usr/bin/env perl

# Get a list of converging / diverging genes

use strict;
use warnings;

# Usage
die "\nUsage: get-gene-orientations.pl <gff-file>\n\n" unless (@ARGV>0);
my $sInput = shift @ARGV;


#########
# START #
#########

# Start processing the gff file
my %hFeatures;
open INPUT, $sInput or die "Can't open $sInput: $!\n";
while (<INPUT>){
   next if (/^\s*$/);
   next if (/^\s*#/);
   next if (lc($_) =~ /^track/);
   s/[\n\r]$//g;
   my ($sChr, $sSource, $sType, $nStart, $nEnd, $nScore, $sStrand, $sFrame, $sID) = split /\t/;
   
   # Check for gff type
   if ($sID =~ /;/){
      $sID = extract_gff2_id('gene_id', $sID);
   }
   
   # Make sure start comes before end
   ($nStart, $nEnd) = sort {$a <=>$b} ($nStart, $nEnd);
   $sChr            = lc($sChr);
   
   # Push results to hash, keep only start and end for ORFs with introns
   if (exists($hFeatures{$sChr}{$sID})){
      $hFeatures{$sChr}{$sID}{'start'}  = $nStart if ($nStart < $hFeatures{$sChr}{$sID}{'start'});
      $hFeatures{$sChr}{$sID}{'end'}    = $nEnd   if ($nEnd   > $hFeatures{$sChr}{$sID}{'end'});
   }
   else{
      $hFeatures{$sChr}{$sID}{'start'}  = $nStart;
      $hFeatures{$sChr}{$sID}{'end'}    = $nEnd;
      $hFeatures{$sChr}{$sID}{'strand'} = $sStrand;
      $hFeatures{$sChr}{$sID}{'type'}   = $sType;
   }
}
close INPUT;


# Start going through each chromosome
foreach my $sChr (keys(%hFeatures)){
   
   # Prepare a new array for processing, needs to have id, strand, start, end
   # Then sort the array descending on everything
   my $rChrFeatures = $hFeatures{$sChr};
   my @aChrFeatures;
   foreach my $sFt (keys(%$rChrFeatures)){
      push @aChrFeatures, [($sFt, $rChrFeatures->{$sFt}{'strand'}, $rChrFeatures->{$sFt}{'start'}, 
                            $rChrFeatures->{$sFt}{'end'}, $rChrFeatures->{$sFt}{'type'})]
   }
   
   # Sort the IDs on the chromosome by position.
   @aChrFeatures = sort {$a->[2] <=> $b->[2] || $a->[3] <=> $b->[3] } @aChrFeatures;
   
   # Now process all entries; when entries overlap, merge them and use the entry with furthest end coordinate to get distance to next feature
   for ( my $i=0 ; $i<@aChrFeatures-1 ; $i++ ){
      my $nCurr = $i;
      for (my $j=$i+1 ; $j<@aChrFeatures ; $j++){
         if ( $aChrFeatures[$nCurr][3] < $aChrFeatures[$j][2]){               # No overlap
            &print_pair($aChrFeatures[$nCurr], $aChrFeatures[$j]);
            $i = $j-1;
            last;
         }
         else{                                                                # Overlap
            $nCurr = $j if ($aChrFeatures[$j][3] > $aChrFeatures[$nCurr][3]);
         }
      }
   }
}

###############
# SUBROUTINES #
###############


sub print_pair{
   my ($rA, $rB) = @_;
   
   my ($sAid, $sAstrand, $nAstart, $nAend, $sAtype) = @$rA;
   my ($sBid, $sBstrand, $nBstart, $nBend, $sBtype) = @$rB;
   my $nDistance = $nBstart-$nAend;
   
   if    ( ($sAstrand eq '+') and ($sBstrand eq '+') ){
      print join("\t", $sAid, '3prime', $sAtype, $sBid, '5prime', $sBtype, 'tandem', $nDistance), "\n";
   }
   elsif ( ($sAstrand eq '+') and ($sBstrand eq '-') ){
      print join("\t", $sAid, '3prime', $sAtype, $sBid, '3prime', $sBtype, '3prime converging', $nDistance), "\n";
   }
   elsif ( ($sAstrand eq '-') and ($sBstrand eq '+') ){
      print join("\t", $sAid, '5prime', $sAtype, $sBid, '5prime', $sBtype, '5prime diverging', $nDistance), "\n";
   }
   elsif ( ($sAstrand eq '-') and ($sBstrand eq '-') ){
      print join("\t", $sAid, '5prime', $sAtype, $sBid, '3prime', $sBtype, 'tandem', $nDistance), "\n";
   }
   else{
      print STDERR "Unknown strand type used for id's '$sAid' and/or '$sBid'\n";
   }
}


sub extract_gff2_id{
   my ($sGffKey, $sGffString)  = @_;
   my @asFields                = split /\s*;\s*/, $sGffString;
   my $sResult                 = '';
   foreach my $sField (@asFields){
      my ($sKey, $sValue) = split /\s"/, $sField;
      $sValue =~ s/"//g;
      $sResult = $sValue if (lc($sKey) eq lc($sGffKey));
   }
   return $sResult;
}
