#!/usr/bin/env perl

# 12.11.2012 11:18:25 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $sTaxID2rank  = "";
my $sName2TaxID  = "";
my $sMegan       = "";
my $sSelectRank  = "phylum";
GetOptions("taxid_mappings:s" => \$sTaxID2rank,
           "name_mappings:s"  => \$sName2TaxID,
           "megan:s"          => \$sMegan,
           "rank:s"           => \$sSelectRank);

# PRINT HELP
$sHelp = 1 unless($sTaxID2rank and $sName2TaxID and $sMegan and $sSelectRank);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
   
   Arguments:
    -t -taxid_mappings <string>
     Name of file with mapping of NCBI taxonomy IDs to rank
     This file can be downloaded from NCBI (nodes.dmp in taxdump.tar.gz)
    -n -name_mappigns <string>  
     Name of file with mapping of taxonomy name to taxonomy ID
     This file must be extracted from the MEGAN data.jar file so that
     the name mappings match the taxonomy IDs exactly.
    -m -megan <string>
     The megan paths export file
    -r -rank <string>
     The rank to extract, e.g. one of phylum, order, family, 
     genus, species
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Fix case
$sSelectRank = lc($sSelectRank);

# Read name to taxonomy ID mappings
my %hName2TaxID;
open TAX, $sName2TaxID or die "Error: can't open '$sName2TaxID': $!\n";
while (<TAX>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my (@asLine) = split /\t/;
   if (@asLine >= 2){
      if (exists $hName2TaxID{$asLine[1]}){
         die "Error: non-unique taxonomy ID '$hName2TaxID{$asLine[1]}' in '$sName2TaxID' on line $.\n";
      }
      else{
         $hName2TaxID{$asLine[1]} = $asLine[0];
      }
   }
   else{
      die "Error: not enough colums in '$sName2TaxID' line $.\n";
   }
}
close TAX;


# Read taxonomy ID to rank mappings
my %hTaxID2rank;
open RANK, $sTaxID2rank or die "Error: can't open '$sTaxID2rank': $!\n";
while (<RANK>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my (@asLine) = split /\t/;
   if (@asLine >= 2){
      if (exists $hTaxID2rank{$asLine[0]}){
         die "Error: non-unique taxonomy ID '$hTaxID2rank{$asLine[1]}' in '$sTaxID2rank' on line $.\n";
      }
      else{
         $hTaxID2rank{$asLine[0]} = lc($asLine[1]);
      }
   }
   else{
      die "Error: not enough colums in '$sTaxID2rank' line $.\n";
   }
}
close RANK;


# Start processing the megan output
open MEGAN, $sMegan or die "Error: can't open '$sMegan': $!\n";
while (<MEGAN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sID, $sEmpty, $sRoot, $nRootScore, $sType, $nTypeScore, $sSuperKingdom, $nSuperKingdomScore, @asPath) = split /\t/;
   
   # Check first for default megan stuff and bacterial reads
   my ($sOutRank, $nOutScore) = ("",0);
   if (lc($sRoot) eq "no hits"){
      $sOutRank  = "unclassified";
      $nOutScore = 100;
   }
   elsif (lc($sRoot) eq "low complexity"){
      $sOutRank  = "lowComplexity";
      $nOutScore = 100;
   }
   elsif (lc($sSuperKingdom) ne "bacteria"){
      $sOutRank  = "noBact";
      $nOutScore = 100;
   }
   else{ # Now process the paths to search for the presence of the rank of interest
      # Convert path array to score hash
      my %hPath;
      for (my $i=0 ; $i<@asPath-2 ; $i+=2){
         die "Error: non-numeric score '$asPath[$i+1]' found on line $. in '$sMegan'\n" unless ($asPath[$i+1] =~ /^\d+$/);
         if (exists $hName2TaxID{$asPath[$i]}){
            my $nTaxID = $hName2TaxID{$asPath[$i]};
            if (exists $hTaxID2rank{$nTaxID}){
               my $sRank = $hTaxID2rank{$nTaxID};
               if (exists $hPath{$sRank}){
                  warn "Warning: found duplicate rank '$sRank' for read '$sID'\n" unless (lc($sRank) eq 'no rank');
               }
               else{
                  $hPath{$sRank}{name}  = $asPath[$i];
                  $hPath{$sRank}{score} = $asPath[$i+1];
               }
            }else{
               warn "Warning: could not find taxID to rank mapping for read '$sID', taxID '$nTaxID' and rank '$asPath[$i]'. Status set to 'unclassified;\n";
            }
         }
         else{
            warn "Warning: could not find name to taxID mapping for read '$sID' and rank '$asPath[$i]'. Status set to 'unclassified'\n";
         }
      }
      
      # Look for rank
      if (exists $hPath{$sSelectRank}){
         $sOutRank  = $hPath{$sSelectRank}{name};
         $nOutScore = $hPath{$sSelectRank}{score};
      }
      else{
         $sOutRank  = "unclassified";
         $nOutScore = 100;
      }
   }
   
   # Print the output to file (after cleaning up the rank name)
   $sOutRank =~ s/\s+\<.*\>//g;
   print join("\t", $sID, $sOutRank, $nOutScore), "\n";
}
close MEGAN;
