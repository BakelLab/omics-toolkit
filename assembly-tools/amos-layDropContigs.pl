#!/usr/bin/env perl

# 10.09.2012 11:25:34 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $sLayoutFile  = "";
my $sIDfile      = "";
GetOptions("help!"   => \$sHelp,
           "lay:s"   => \$sLayoutFile,
           "id:s"    => \$sIDfile);

# PRINT HELP
$sHelp = 1 unless($sLayoutFile and (@ARGV or $sIDfile));
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -l <lay file> [-i <id file>] [id1 .. idN]
   
   Removes LAY entries for supplied contig identifiers so that
   they are not used to build amos banks. The sequences associated
   with the contigs (i.e. the RED entries) are kept as 'orphaned'
   entries.
   
   Options:
    -l --lay <string>
      Layout file. LAY entries with matching eids will be dropped.
    -i --ids <string>
      Optional file with ids to drop.
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Read IDS from file
my %hDropIDs;
if ($sIDfile){
   open IN, "<$sIDfile" or die "Error: can't read ID file\n";
   while (<IN>){
      next if (/^\s*$/);
      next if (/^ *#/);
      s/[\n\r]+$//;
      my (@asIDs) = split /\s/;
      foreach my $sID (@asIDs){
         $hDropIDs{$sID} = "";
      }
   }
   close IN;
}

# Read IDs from command line
foreach my $sID (@ARGV){
   $hDropIDs{$sID} = "";
}

# Process the layout file
my $nDropCount = 0;
my ($sBuffer, $sEid) = ("","");
open IN, $sLayoutFile or die "Error: can't open '$sLayoutFile': $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   if(/^{LAY/){
      if(exists $hDropIDs{$sEid}){
         $nDropCount++;
      }
      else{
         print $sBuffer;
      }
      $sBuffer = $_;
      next;
   }
   if(/^{RED/){
      if(exists $hDropIDs{$sEid}){
         $nDropCount++;
      }
      else{
         print $sBuffer;
      }
      $sBuffer = "";
   }
   if($sBuffer){
      $sBuffer .= $_;
   }
   else{
      print;
   }
   $sEid = $1 if (/^eid:(.*)$/);
}
close IN;

# Status message
my $sMessage = $nDropCount == 1 ? "Dropped 1 contig layout\n" : "Dropped $nDropCount contigs\n";
warn $sMessage;
