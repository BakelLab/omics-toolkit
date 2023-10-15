#!/usr/bin/env perl

# 08.01.2014 16:18:26 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp         = 0;
my $sRefFeature   = 'exon';
my $sAttributes   = "transcript_id";
my $sJoinChar     = "|";
GetOptions("help!"       => \$sHelp,
           "feature:s"   => \$sRefFeature,
           "attribute:s" => \$sAttributes,
           "join:s"      => \$sJoinChar);

# GLOBALS
my %hMultiCopyTags = ('tag'=>'','Dbxref'=>'','gbkey'=>'','Name'=>'','product'=>'','Note'=>'','tss_id'=>'', 'ont'=>'', 'db_xref'=>'', 'inference'=>'');

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName <gtf-file>
   
   Convert a gtf file with an attribute list to a basic gff file
   that has a single attribute linking lines corresponding to one
   feature/item.
   
   Arguments:
    -f --feature <string>
      Comma-separated list of features to extract mappings from.
      Default: $sRefFeature
    -a --attribute <string>
      The attribute tag(s) to use as group ID in the gff file. 
      Multiple attributes can be combined into a compound ID by 
      specifying a comma-separated list of attribute tags.
      Default: $sAttributes
    -j --join <string>
      Character to join the attribute tags with. Default: $sJoinChar
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Parse attribute list
my @asAttributes = split /\,/, $sAttributes;

# Parse feature list
my %hFeatures;
for my $sFeature (split /\,/, lc($sRefFeature)){
   $hFeatures{$sFeature}++;
}

# Process gtf
open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]': $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sSeqName, $sSource, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sGroup) = split /\t/, $_, -1;
   my $rhAnnots = gtf_annots_to_hash($sGroup, \%hMultiCopyTags);
   
   # Look for attributes
   my $sGffID    = '';
   my $flMissing = 1;
   foreach my $sAttrib (@asAttributes){
      if(exists $rhAnnots->{$sAttrib}){
         my $sCompoundAttrib = exists($hMultiCopyTags{$sAttrib}) ? join($sJoinChar, @{$rhAnnots->{$sAttrib}}) : $rhAnnots->{$sAttrib};
         $sGffID .= $sGffID ? $sJoinChar . $sCompoundAttrib : $sCompoundAttrib;
         $flMissing = 0;
      }
   }
   
   # Figure out if this is the feature type we want
   my $flFeatureMatch = 0;
   if ( ($sRefFeature eq "") or (exists $hFeatures{lc($sFeature)}) ){
      $flFeatureMatch = 1;
   }
   
   # Print the gff line
   if ($flFeatureMatch){
      if ( !$flMissing ){
         print join("\t", $sSeqName, $sSource, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sGffID), "\n";
      }
      else{
         warn "Warning: skipping gtf line $.: missing attributes\n";
      }
   }
}
close IN;

#################
## SUBROUTINES ##
#################

# gtf_annots_to_hash
#
# Parse gtf annotations and return key-value pairs
sub gtf_annots_to_hash {
   my ($sAnnots, $rMultiCopy) = @_;
   my %hReturn;
   
   my @asPairs = split /" *; */, $sAnnots;
   foreach my $sPair (@asPairs){
      my ($sKey, $sVal) = split / "/, $sPair;
      if ( exists $rMultiCopy->{$sKey} ){
         push @{$hReturn{$sKey}}, $sVal;
      }
      else{
         die "Error: Duplicate key entry '$sKey' found with values '$hReturn{$sKey}' and '$sVal'\n" if (exists $hReturn{$sKey});
         $hReturn{$sKey} = $sVal;
      }
   }
   return \%hReturn;
}
