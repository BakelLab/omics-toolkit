#!/usr/bin/env perl

# 08.01.2014 15:09:37 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GLOBALS
my %hMultiCopyTags = ('tag'=>'','Dbxref'=>'','gbkey'=>'','Name'=>'','product'=>'','Note'=>'','tss_id'=>'', 'ont'=>'', 'db_xref'=>'', 'gene_synonym'=>'');

# GET PARAMETERS
my $sHelp        = 0;
my $flInvert     = 0;
my $sMatchFile   = "";
my $sAttribute   = 'transcript_id';
my $sAddReplace  = 'gene_id';
my $sDefault     = '';
GetOptions("help!"        => \$sHelp,
           "match:s"      => \$sMatchFile,
           "attribute:s"  => \$sAttribute,
           "replace:s"    => \$sAddReplace,
           "default:s"    => \$sDefault);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0 and $sMatchFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [-a] -m <match-file> <gtf-file>
   
   Filter a gtf file based on attribute matches, e.g.
   gene_id or gene_type
   
   Arguments:
    -m --match <string>
      Attribute add/replace table in tab delimited format.
      First column should contain attribute values to match on.
      The second column should contain the values of the attribute
      to add or replace.
    -a --attribute <string>
      Name of attribute to match on
      Default: $sAttribute
    -r --replace <string>
      Name of attribute to add/replace
      Default: $sAddReplace
    -d --default <string>
      Default value to use when the attribute tag does not exist
      and there is no match in the add/replace table.
      This argument ensures that the add/replace tag is always
      set to a default value.
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# We don't support the tag attribute as it can occur more than once
die "Error: 'tag' attribute is unsupported\n" if ($sAttribute eq 'tag');

# Read ID list
my %hLUT;
open MATCH, $sMatchFile or die "Error: can't open '$sMatchFile': $!\n";
while (<MATCH>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my (@asPair) = split /\s+/;
   if (@asPair == 2){
      my ($sKey, $sVal) = @asPair;
      if (exists $hLUT{$sKey}){
         die "Error: duplicate entry '$sKey' in match file, exiting\n";
      }
      else{
         $hLUT{$sKey} = $sVal;
      }
   }
   else{
      warn("Skipping '$_': missing valid key/value pairs\n");
   }
}
close MATCH;


# Match gtf attribute
my ($nMatchCount, $nNoMatchCount, $nNotFoundCount) = (0, 0, 0);
open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]': $!\n";
while (<IN>){
   next if (/^\s*$/);
   if (/^ *#/){
      print;
      next;
   };
   s/[\n\r]+$//;
   my ($sSeqName, $sSource, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sGroup) = split /\t/, $_, -1;
   my ($rhAnnots, $raKeyOrder) = gtf_annots_to_hash($sGroup, \%hMultiCopyTags);
   if (exists $rhAnnots->{$sAttribute}){
      my $sKey = $rhAnnots->{$sAttribute};
      if (exists $hLUT{$sKey} ){
         push @{$raKeyOrder}, $sAddReplace unless ( exists $rhAnnots->{$sAddReplace} ); # Append the new key at the end of the gtf string if it doesn't already exist!
         $rhAnnots->{$sAddReplace} = $hLUT{$sKey};
         print join ("\t", $sSeqName, $sSource, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $nFrame, hash_to_gtf_annots($rhAnnots, $raKeyOrder, \%hMultiCopyTags) ), "\n";
         $nMatchCount++;
      }
      else{
         if ( exists $rhAnnots->{$sAddReplace} ){
            warn("Warning: did not replace the $sAddReplace tag for $sAttribute '$sKey' because it was not found in the match file\n");
            print "$_\n" ;
         }
         else {
            if ($sDefault){
               push @{$raKeyOrder}, $sAddReplace;
               $rhAnnots->{$sAddReplace} = $sDefault;
               warn("Warning: set the $sAddReplace tag to default value '$sDefault' for $sAttribute '$sKey' because it was not found in the match file\n");
               print join ("\t", $sSeqName, $sSource, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $nFrame, hash_to_gtf_annots($rhAnnots, $raKeyOrder, \%hMultiCopyTags) ), "\n";
            }
            else{
               warn("Warning: did not set the $sAddReplace tag for $sAttribute '$sKey' because it was not found in the match file\n");
               print "$_\n" ;
            }
         }
         $nNoMatchCount++;
      }
   }
   else{
      if ( exists $rhAnnots->{$sAddReplace} ){
         warn("Warning: did not replace the $sAddReplace tag for gtf line $. because it did not contain the $sAttribute tag\n");
         print "$_\n" ;
      }
      else {
         if ($sDefault){
            push @{$raKeyOrder}, $sAddReplace;
            $rhAnnots->{$sAddReplace} = $sDefault;
            warn("Warning: set the $sAddReplace tag to default value '$sDefault' for gtf line $. because it did not contain the $sAttribute tag\n");
            print join ("\t", $sSeqName, $sSource, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $nFrame, hash_to_gtf_annots($rhAnnots, $raKeyOrder, \%hMultiCopyTags) ), "\n";
         }
         else{
            warn("Warning: did not set the $sAddReplace tag for gtf line $. because it did not contain the $sAttribute tag\n");
            print "$_\n" ;
         }
      }
      $nNotFoundCount++;
   }
}
close IN;

# Report replacement stats
warn "\n";
warn("Added or replaced the '$sAddReplace' tag in $nMatchCount GTF lines containing the '$sAttribute' tag\n") if ($nMatchCount);
if ($sDefault){
   warn("Set the '$sAddReplace' tag to default value '$sDefault' in $nNoMatchCount GTF lines because the '$sAttribute' tag value was not present in the match file\n") if $nNoMatchCount;
   warn("Set the '$sAddReplace' tag to default value '$sDefault' in $nNotFoundCount GTF lines because they did not contain the '$sAttribute' tag\n") if $nNotFoundCount;
}
else{
   warn("$nNoMatchCount GTF lines were unchanged because the '$sAttribute' tag value was not present in the match file\n") if ($nNoMatchCount);
   warn("$nNotFoundCount GTF lines were unchanged because they did not contain the '$sAttribute' tag\n") if ($nNotFoundCount);
}

#################
## SUBROUTINES ##
#################

# gtf_annots_to_hash(gtf_string)
#
# Parse gtf annotations and return key-value pairs
sub gtf_annots_to_hash {
   my ($sAnnots, $rMultiCopy) = @_;
   my %hReturn;
   
   $sAnnots =~ s/;+ *$//;
   $sAnnots =~ s/"//g;
   my @asKeyOrder;
   my @asPairs = split / *; */, $sAnnots;
   foreach my $sPair (@asPairs){
      my ($sKey, @asVal) = split / /, $sPair;
      my $sVal = join(" ", @asVal);
      if (exists $rMultiCopy->{$sKey}){
         push @asKeyOrder, $sKey unless (exists $hReturn{$sKey}); 
         push @{$hReturn{$sKey}}, $sVal;
      }
      else{
         die "Error: Duplicate key entry '$sKey' found with values '$hReturn{$sKey}' and '$sVal'\n" if (exists $hReturn{$sKey});
         $hReturn{$sKey} = $sVal;
         push @asKeyOrder, $sKey;
      }
   }
   return (\%hReturn, \@asKeyOrder);
}

# hash_to_gtf_annots(gtf_hash, key_order_array)
#
# Convert hash of gtf tag-value pairs back into a gtf string
sub hash_to_gtf_annots {
   my ($rhAnnots, $raKeyOrder, $rMultiCopy) = @_;

   my @asReturn;
   foreach my $sKey (@$raKeyOrder){
      if (exists $rMultiCopy->{$sKey}){
         foreach my $sVal (@{$rhAnnots->{$sKey}}){
            push @asReturn, join("", $sKey, ' "', $sVal, '";') unless (($sKey eq "") or ($sVal eq ""));
         }
      }
      elsif ($sKey eq "exon_number"){
         push @asReturn, join("", $sKey, ' ', $rhAnnots->{$sKey}, ';') unless (($sKey eq "") or ($rhAnnots->{$sKey} eq ""));
      }
      else{
         push @asReturn, join("", $sKey, ' "', $rhAnnots->{$sKey}, '";') unless (($sKey eq "") or ($rhAnnots->{$sKey} eq ""));
      }
   }
   return join(" ", @asReturn);
}
