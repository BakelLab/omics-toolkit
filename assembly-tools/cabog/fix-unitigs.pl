#!/usr/bin/perl

# 17.01.2011 09:44:34 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GENERAL
$ENV{TIGSTORE} = "tigStore";  # Path to tigstore binary

# GET PARAMETERS
my $sHelp        = 0;
my $nSegment     = 0;
my $sGkpStore    = "";
my $sTigStore    = "";
GetOptions("help!"       => \$sHelp,
           "up:n"        => \$nSegment,
           "gkp:s"       => \$sGkpStore,
           "tig:s"       => \$sTigStore);
my $sInput = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sInput and $nSegment and $sGkpStore and $sTigStore);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -g <gkpStore> -t <tigStore> -up <segment> <.err file>
   
   Fixes errors in unitig consensus sequence building by
   moving unalignable fragments to a new unitig.
    
   Options: 
    -g <string>
      Path to gatekeeper store
    -t <string>
      Path to tigstore
    -up <integer>
      Segment ID (must match error file segment ID)
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check arguments
die "Error: could not find gatekeepr store\n" unless(-e $sGkpStore);
die "Error: could not find unitig store\n" unless(-e $sTigStore);

# Get a hash with failed fragments
my %hFragmentErrors;
open IN, $sInput or die "Error: can't open '$sInput': $!\n";
while (<IN>){
   next if /^\s*$/;
   if (/Unitig (\d+) FAILED.* Could not align fragment (\d+)/){
      $hFragmentErrors{$1}{$2} = '';
   }
}
close IN;

# Print message with the number of unitig errors found
my $nErrCount = scalar(keys %hFragmentErrors);
my $sErrText  = $nErrCount==1 ? "error" : "errors";
warn(join('', "Found $nErrCount unitig $sErrText\n"));

# Process each unitig
foreach my $nUnitig (keys %hFragmentErrors){
   system("$ENV{TIGSTORE} -g $sGkpStore -t $sTigStore 1 -up $nSegment -d layout -u $nUnitig > unitig$nUnitig") == 0 
      or die "Error: could not export unitig '$nUnitig' from the tigstore: $!";
   
   # Process the unitig exported
   my @asHeader;
   my %hHeader;
   my @asKeepFrags;
   my @asMoveFrags;
   open UNITIG, "unitig$nUnitig" or die "Error: can't process unitig $nUnitig\n";
   while (<UNITIG>){
      next if (/^\s*$/);
      chomp;
      my @asLine = split /\s+/;
      if    ($asLine[0] eq "unitig")  { $hHeader{$asLine[0]} = $asLine[1]; push @asHeader, $asLine[0];}
      elsif ($asLine[0] eq "len")     { $hHeader{$asLine[0]} = $asLine[1]; push @asHeader, $asLine[0];}
      elsif ($asLine[0] eq "cns")     { $hHeader{$asLine[0]} = '';         push @asHeader, $asLine[0];}
      elsif ($asLine[0] eq "qlt")     { $hHeader{$asLine[0]} = '';         push @asHeader, $asLine[0];}
      elsif ($asLine[0] =~ "^data\.") { $hHeader{$asLine[0]} = $asLine[1]; push @asHeader, $asLine[0];}
      elsif ($asLine[0] eq "FRG") { 
         if(exists $hFragmentErrors{$nUnitig}{$asLine[4]}){ push @asMoveFrags, $_ }
         else{ push @asKeepFrags, $_ }
      }
      else{
         die "Error: Uknown field in unitig: $asLine[0]";
      }
   }
   close UNITIG;
   die "Error: Mismatch" unless (scalar(@asMoveFrags) == scalar(keys(%{$hFragmentErrors{$nUnitig}})));
   die "Error: missing the data.num_frags field in the unitig header" unless (exists $hHeader{'data.num_frags'});

   # Print the two new unitigs
   open FIX, ">unitig$nUnitig.fix" or die "Error: write fix file for unitig $nUnitig\n";
   foreach my $sKey (@asHeader){
      if ($sKey eq "data.num_frags"){
         print FIX join(" ", $sKey, scalar(@asKeepFrags) ), "\n";
      }
      else{
         print FIX join(" ", $sKey, $hHeader{$sKey}), "\n";
      }
   }
   foreach my $sFrg (@asKeepFrags){
      print FIX $sFrg, "\n";
   }
   foreach my $sKey (@asHeader){
      if ($sKey eq "data.num_frags"){
         print FIX join(" ", $sKey, scalar(@asMoveFrags) ), "\n";
      }
      elsif ($sKey eq "unitig"){
         print FIX join(" ", $sKey, -1), "\n";
      }
      else{
         print FIX join(" ", $sKey, $hHeader{$sKey}), "\n";
      }
   }
   foreach my $sFrg (@asMoveFrags){
      print FIX $sFrg, "\n";
   }
   close FIX;

   # Finally re-insert the updated unitig
   system("$ENV{TIGSTORE} -g $sGkpStore -t $sTigStore 1 -up $nSegment -R unitig$nUnitig.fix") == 0 
      or die "Error: could not re-insert unitig '$nUnitig' into the tigstore: $!";
}

# Rerun consensus after all the unitigs have been updated
system("sh consensus.sh $nSegment") == 0
   or die "Error: could not rerun consensus: $!";
