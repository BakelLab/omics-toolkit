#!/usr/bin/perl

# 24.11.2014 15:33:11 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;

# GET PARAMETERS
my $sHelp        = 0;
GetOptions("help!"   => \$sHelp);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

my $seqio_object = Bio::SeqIO->new(-file => $ARGV[0] );
while (my $seq_object = $seqio_object->next_seq){
   my ($sAccession, $sStrain, $sOrganism, $sSeq, $sProtName) = ("Unk","Unk","Unk","Unk","Unk");
   $sAccession = $seq_object->accession();

   # Cycle through sequence features
   for my $feat_object ($seq_object->get_SeqFeatures) {
      if (lc($feat_object->primary_tag) eq "source"){
         my %feat_tags = map {$_ => 1} ($feat_object->get_all_tags);
         ($sStrain)   =  ($feat_object->get_tag_values('strain')) if (exists $feat_tags{'strain'});
         ($sOrganism) =  ($feat_object->get_tag_values('organism')) if (exists $feat_tags{'organism'});
         $sOrganism   = lc($sOrganism);
         $sOrganism   =~ s/peptoclostridium/C/;
         $sOrganism   =~ s/clostridium/C/;
         $sOrganism   =~ s/difficile/diff/;
      }
      if (lc($feat_object->primary_tag) eq "protein"){
         my %feat_tags = map {$_ => 1} ($feat_object->get_all_tags);
         ($sProtName)   =  ($feat_object->get_tag_values('name')) if (exists $feat_tags{'name'});
         ($sProtName)   =  ($feat_object->get_tag_values('product')) if (exists $feat_tags{'product'});
      }
   }
   my $sFastaHeader = ">${sAccession}|$sProtName|${sOrganism}_${sStrain}";
   $sFastaHeader =~ s/ /_/g;
   my $sFastaSeq = $seq_object->seq;
   $sFastaSeq =~ s/.{100}/$&\n/sg;
   $sFastaSeq =~ s/\n+$//;
   print "$sFastaHeader\n$sFastaSeq\n";
}
