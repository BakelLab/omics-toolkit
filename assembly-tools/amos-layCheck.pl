#!/usr/bin/env perl

# 10.09.2012 11:25:34 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
GetOptions("help!"   => \$sHelp);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName <.lay file>
   
   Checks an amos layout file (such as those produced by the
   cabog pacBioToCA pipeline) for a variety of errors.
   
   Arguments: 
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

my @asLayouts;
my %hReadInfo;

my $sLastElement = "";
my $nContigID    = 0;
my ($sEid, $sIid,$sOff, $sSrc, $nCls, $nCle, $nSeq, $nQlt) = ("","",""."",0,0,0,0);
open IN, $ARGV[0] or die "Error: can't open input\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   if(/LAY/ or /TLE/ or /RED/){
      if($sLastElement eq '{LAY'){
         $nContigID++;
      }
      elsif($sLastElement eq '{TLE'){
         push @asLayouts, [($nContigID, $sEid, $sIid, $nCls, $nCle, $sOff, $sSrc)];
      }
      elsif($sLastElement eq '{RED'){
         $hReadInfo{$sEid} = [($nCls, $nCle,$nSeq,$nQlt)];
      }
      else{
      }
      $sLastElement = $_;
   }
   $sEid = $1 if (/^eid:(.*$)/);
   $sIid = $1 if (/^iid:(.*$)/);
   $sOff = $1 if (/^off:(.*$)/);
   $sSrc = $1 if (/^src:(.*$)/);
   ($nCls, $nCle) = split /,/, $1 if (/^clr:(.*$)/);
   if (/^seq:$/){
      $_ = <IN>;
      $nSeq = length($_);
   }
   if (/^qlt:$/){
      $_ = <IN>;
      $nQlt = length($_);
   }
}
close IN;

foreach my $rLayout (@asLayouts){
   my $sID = pop @$rLayout;
   if(exists $hReadInfo{$sID}){
      my ($nContigID, $sEid, $sIid, $nTLEcls, $nTLEcle, $sOff) = @$rLayout;
      my ($nREDcls, $nREDcle,$nSeq,$nQlt) = @{$hReadInfo{$sID}};
      # Check for qual and seq length diff
      print join("\t", "Seq/Qual diff:", @$rLayout, $sID, @{$hReadInfo{$sID}}), "\n" unless($nSeq == $nQlt);
      
      # Check for offset length diff
      my $nTLEdiff = abs($nTLEcls - $nTLEcle);
      my $nREDdiff = abs($nREDcls - $nREDcle);
      print join("\t", "Offset diff:", @$rLayout, $sID, @{$hReadInfo{$sID}}), "\n" unless($nTLEdiff == $nREDdiff);
      
      # Check for range errors
      print join("\t", "Range error:", @$rLayout, $sID, @{$hReadInfo{$sID}}), "\n" if ($nTLEcls>=$nSeq or $nTLEcle>=$nSeq or $nREDcls>=$nSeq or $nREDcle>=$nSeq);
   }
   else{
      print join("\t", "Missing RED:", @$rLayout, $sID), "\n";
   }
}
