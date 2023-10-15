#!/usr/bin/env perl

# 08.01.2014 15:09:37 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $flInvert     = 0;
my $sMatchFile   = "";
my $sAttribute   = 'transcript_id';
GetOptions("help!"        => \$sHelp,
           "match:s"      => \$sMatchFile,
           "attribute:s"  => \$sAttribute,
           "v!"           => \$flInvert);
           
# GLOBALS
my %hMultiCopyTags = ('tag'=>'','Dbxref'=>'','gbkey'=>'','Name'=>'','product'=>'','Note'=>'','tss_id'=>'', 'ont'=>'', 'db_xref'=>'', 'gene_synonym'=>'');

# PRINT HELP
$sHelp = 1 unless(@ARGV>0 and $sMatchFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [-a] -m <match-file> <gtf-file>
   
   Filter a gtf file based on attribute matches, e.g.
   gene_id or gene_type
   
   Arguments:
    -a --attribute <string>
      Attribute name to filter on
      Default: $sAttribute
    -m --match <string>
      Name of file containing IDs to filter on
    -v
      Invert match (i.e. select all entries NOT matching
      the attribute criteria)
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Read ID list
my %hMatchIDs;
open MATCH, $sMatchFile or die "Error: can't open '$sMatchFile': $!\n";
while (<MATCH>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my @asMatch = split /\s+/;
   foreach my $sMatch (@asMatch) {$hMatchIDs{$sMatch}++}
}
close MATCH;

# Match gtf attribute
my $nMatchCount = 0;
open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]': $!\n";
while (<IN>){
   next if (/^\s*$/);
   if (/^ *#/){
      print;
      next;
   };
   s/[\n\r]+$//;
   my ($sSeqName, $sSource, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sGroup) = split /\t/, $_, -1;
   my $rhAnnots = gtf_annots_to_hash($sGroup, \%hMultiCopyTags);
   my $flMatch = 0;
   if (exists $rhAnnots->{$sAttribute}){
      if (exists $hMatchIDs{$rhAnnots->{$sAttribute}}){
         $flMatch = 1;
      }
   }
   if ($flInvert){
      unless ($flMatch){
         print "$_\n" ;
         $nMatchCount++;
      }
   }
   else{
      if ($flMatch){
         print "$_\n";
         $nMatchCount++;
      }
   }
}
close IN;

warn("Extracted $nMatchCount lines\n");


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
