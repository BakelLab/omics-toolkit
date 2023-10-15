#!/usr/bin/env perl

# 08.01.2014 15:09:37 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $sRefFeature  = 'exon';
my $sAttributes  = 'gene_id,gene_name';
GetOptions("help!"        => \$sHelp,
           "feature:s"    => \$sRefFeature,
           "attributes:s" => \$sAttributes);

# GLOBALS
my %hMultiCopyTags = ('tag'=>'','Dbxref'=>'','gbkey'=>'','Name'=>'','product'=>'','Note'=>'','tss_id'=>'', 'ont'=>'', 'db_xref'=>'', 'gene_synonym'=>'');

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [-sa] <gtf-file>
   
   Count how many times a particular attribute value or combination
   of attribute values occurs in a gtf file. Can also be used to
   extract between-attribute mappings, e.g. gene_id to gene_name.
   
   Arguments:
    -f --feature <string>
      Source feature to extract mappings from.
      Default: $sRefFeature
    -a --attributes <string>
      Comma-separated list of one or more attributes to extract.
      Default: $sAttributes
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Parse attributes
my @asAttributes = split /,/, $sAttributes;

# Collect attribute combination counts
my %hOut;
open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]': $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sSeqName, $sSource, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sGroup) = split /\t/, $_, -1;
   if ($sFeature eq $sRefFeature){
      my $rhAnnots = gtf_annots_to_hash($sGroup, \%hMultiCopyTags);
      my @asKeys;
      foreach my $sAttrib (@asAttributes){
         my $sKey = '';
         if ( exists $rhAnnots->{$sAttrib} ){
            $sKey = exists($hMultiCopyTags{$sAttrib}) ? join("|", @{$rhAnnots->{$sAttrib}}) : $rhAnnots->{$sAttrib};
         }
         push @asKeys, $sKey;
      }
      $hOut{join("\t", @asKeys)}{featurecount}++;
      
      if (exists $rhAnnots->{'transcript_id'}){
         $hOut{join("\t", @asKeys)}{transcriptcount}{$rhAnnots->{'transcript_id'}}++;
      }
      if (exists $rhAnnots->{'gene_id'}){
         $hOut{join("\t", @asKeys)}{genecount}{$rhAnnots->{'gene_id'}}++;
      }
   }
}
close IN;

# Write output
print "#", join("\t", @asAttributes, "${sRefFeature}count", "transcriptcount", "genecount"), "\n";
foreach my $sKey (sort keys %hOut){
   my $nTranscriptCount = exists($hOut{$sKey}{transcriptcount}) ? scalar(keys %{$hOut{$sKey}{transcriptcount}}) : 0;
   my $nGeneCount       = exists($hOut{$sKey}{genecount})       ? scalar(keys %{$hOut{$sKey}{genecount}}) : 0;
   print join("\t", $sKey, $hOut{$sKey}{featurecount}, $nTranscriptCount, $nGeneCount), "\n";
}

#################
## SUBROUTINES ##
#################

# gtf_annots_to_hash
#
# Parse gtf annotations and return key-value pairs
sub gtf_annots_to_hash {
   my ($sAnnots, $rMultiCopy) = @_;
   my %hReturn;
   
   my @asPairs = split / *; */, $sAnnots;
   foreach my $sPair (@asPairs){
      my ($sKey, $sVal) = split(/ /, $sPair, 2);
      $sVal =~ s/"//g;
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
