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

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [-sa] <gff-file>
   
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
      my $rhAnnots = gtf_annots_to_hash($sGroup);
      my @asKeys;
      foreach my $sAttrib (@asAttributes){
         my $sKey = exists $rhAnnots->{$sAttrib} ? $rhAnnots->{$sAttrib} : '';
         push @asKeys, $sKey;
      }
      $hOut{join("\t", @asKeys)}++;
   }
}
close IN;

# Write output
print "#", join("\t", @asAttributes, "count"), "\n";
foreach my $sKey (sort keys %hOut){
   print join("\t", $sKey, $hOut{$sKey}), "\n";
}

#################
## SUBROUTINES ##
#################

# gtf_annots_to_hash
#
# Parse gtf annotations and return key-value pairs
sub gtf_annots_to_hash {
   my ($sAnnots) = @_;
   my %hReturn;
   
   my @asPairs = split / *; */, $sAnnots;
   foreach my $sPair (@asPairs){
      my ($sKey, $sVal) = split /=/, $sPair;
      $sVal =~ s/"//;
      if ($sKey eq 'tag'){
         push @{$hReturn{$sKey}}, $sVal;
      }
      else{
         die "Error: Duplicate key entry '$sKey' found with values '$hReturn{$sKey}' and '$sVal'\n" if (exists $hReturn{$sKey});
         $hReturn{$sKey} = $sVal;
      }
   }
   return \%hReturn;
}
