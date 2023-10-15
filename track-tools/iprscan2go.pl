#!/usr/bin/env perl

# 15.10.2012 13:56:22 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

# GET PARAMETERS
my $sHelp        = 0;
my $sIprsFile    = "";
my $sGOtermsFile = "";
my $sGeneIdFile  = "";
my $sDBDate      = strftime "%Y%m%d", localtime;
my $sDBTaxon     = '79327';
my $sDB          = "SGD";
my $sDBref       = "PMID:0000001";
my $sDBEvi       = 'IEA';
GetOptions("help!"    => \$sHelp,
           "scan:s"   => \$sIprsFile,
           "term:s"   => \$sGOtermsFile,
           "ids:s"    => \$sGeneIdFile,
           "adate:s"  => \$sDBDate,
           "xtaxon:s" => \$sDBTaxon,
           "db:s"     => \$sDB,
           "ref:s"    => \$sDBref,
           "evi:s"    => \$sDBEvi);

# PRINT HELP
$sHelp = 1 unless($sIprsFile and $sGOtermsFile and $sGeneIdFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -s <iprs-file> -t <GO-terms-file> -i <geneID-file>
   
   Script to convert an interproscan file to a gaf formatted
   GO reference file to use for GO term enrichment analysis.
   
   Options:
    -s <string> 
      InterproScan output file (raw format + Interpro + GO lookup)
    -t <string>
      GO terms file, download from http://www.geneontology.org/GO.downloads.files.shtml
    -i <string>
      File with all gene identifiers in the reference genome
      Gene IDs without interproscan hits will be listed as "protein of unknown function"
    -a <numeric>
      Date in yyyymmdd format. Default: current date
    -x <numeric>
      NCBI taxonomy ID for organism. Default: $sDBTaxon
    -r <string>
      Annotation source reference. Default: $sDBref
    -e <string>
      GO evidence code. Default: $sDBEvi
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Read GO terms
my %hGOterms;
open GO, $sGOtermsFile or die "Error: can't open $sGOtermsFile: $!\n";
while (<GO>){
   next if (/^\s*$/);
   next if (/^ *#/);
   next if (/^\!/);
   s/[\n\r]+$//;
   my ($sGOid, $sAltGOid, $sName, $sOntology, $sObsolete) = split /\t/, $_, -1;
   $hGOterms{$sGOid}{alt}      = $sAltGOid;
   $hGOterms{$sGOid}{name}     = $sName;
   $hGOterms{$sGOid}{ontology} = $sOntology;
   $hGOterms{$sGOid}{obsolete} = 1 if ($sObsolete);
}
close GO;

# Read InterproScan hits
my %hIprsHits;
open IN, $sIprsFile or die "Error: can't open $sIprsFile: $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sSeqID, $sCheck, $nNo, $sSrcDb, $sSrcID,$sSrcDescr, $nSrcStart, $nSrcEnd, $nSrcE, $sSrcCode, $sDate, $sIPR, $sIPRDesc, $sGOfield) = split /\t/;
   my @asGO = ($sGOfield =~ /(GO:\d+)/g);
   foreach my $sGO (@asGO){$hIprsHits{$sSeqID}{$sGO}++}
}
close IN;

# Read geneIDs and process gene ontology stuff
my $sDBObjectType = 'gene';
my $sDBAssignedBy = $sDB ;
open ID, $sGeneIdFile or die "Error: can't open $sGeneIdFile: $!\n";
print "!gaf-version: 2.0\n";
while (<ID>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my %hOntology;
   if (exists $hIprsHits{$_}){
      foreach my $sGOid (keys %{$hIprsHits{$_}}){
         print join("\t", $sDB, $_, $_, "", $sGOid, $sDBref, $sDBEvi,"", $hGOterms{$sGOid}{ontology},"","",$sDBObjectType,"taxon:$sDBTaxon", $sDBDate, $sDBAssignedBy), "\n";
         $hOntology{$hGOterms{$sGOid}{ontology}}++;
      }
   }
   print join("\t", $sDB, $_, $_, "", "GO:0005575", $sDBref, "ND", "", "C","","",$sDBObjectType,"taxon:$sDBTaxon", $sDBDate, $sDBAssignedBy), "\n" unless(exists $hOntology{'C'});
   print join("\t", $sDB, $_, $_, "", "GO:0003674", $sDBref, "ND", "", "F","","",$sDBObjectType,"taxon:$sDBTaxon", $sDBDate, $sDBAssignedBy), "\n" unless(exists $hOntology{'F'});
   print join("\t", $sDB, $_, $_, "", "GO:0008150", $sDBref, "ND", "", "P","","",$sDBObjectType,"taxon:$sDBTaxon", $sDBDate, $sDBAssignedBy), "\n" unless(exists $hOntology{'P'});
}
close ID;


