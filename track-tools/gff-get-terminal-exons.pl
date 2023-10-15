#!/usr/bin/env perl

# 16.08.2018 16:34:57 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $flSkipMono   = 0;
GetOptions("mono!"   => \$flSkipMono,
           "help!"   => \$sHelp);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [-m] gff-file
    
   Arguments:
    -m
      Skip mono-exon genes
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

my @aOrder;
my %hExons;

open IN, $ARGV[0] or die "Error: can't open $ARGV[0]: $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sChr, $sSrc, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sName) = split(/\t/, $_, -1);
   if (exists $hExons{$sName}){
      if ( ($hExons{$sName}{chr} eq $sChr) and ($hExons{$sName}{strand} eq $sStrand)){
         $hExons{$sName}{distinct_exons}{"$sChr|$nStart|$nEnd|$sStrand"}++;
         if ($sStrand eq "+"){
            if ( ($nStart < $hExons{$sName}{first}{start}) or ( ($nStart == $hExons{$sName}{first}{start}) and ($nEnd > $hExons{$sName}{first}{end}) ) ){
               $hExons{$sName}{first}{start} = $nStart;
               $hExons{$sName}{first}{end}   = $nEnd;
               $hExons{$sName}{first}{score} = $nScore;
               $hExons{$sName}{first}{frame} = $nFrame;
            }
            if ( ($nEnd > $hExons{$sName}{last}{end}) or ( ($nEnd == $hExons{$sName}{last}{end}) and ($nStart < $hExons{$sName}{last}{start}) ) ){
               $hExons{$sName}{last}{start} = $nStart;
               $hExons{$sName}{last}{end}   = $nEnd;
               $hExons{$sName}{last}{score} = $nScore;
               $hExons{$sName}{last}{frame} = $nFrame;
            }
         }
         else{         
            if ( ($nEnd > $hExons{$sName}{first}{end}) or ( ($nEnd == $hExons{$sName}{first}{end}) and ($nStart < $hExons{$sName}{first}{start}) ) ){
               $hExons{$sName}{first}{start} = $nStart;
               $hExons{$sName}{first}{end}   = $nEnd;
               $hExons{$sName}{first}{score} = $nScore;
               $hExons{$sName}{first}{frame} = $nFrame;
            }
            if ( ($nStart < $hExons{$sName}{last}{start}) or ( ($nStart == $hExons{$sName}{last}{start}) and ($nEnd > $hExons{$sName}{last}{end}) ) ){
               $hExons{$sName}{last}{start} = $nStart;
               $hExons{$sName}{last}{end}   = $nEnd;
               $hExons{$sName}{last}{score} = $nScore;
               $hExons{$sName}{last}{frame} = $nFrame;
            }
         }   
      }
      else{
         warn("Feature '$sName' was annotated for more than one chromosome or strand. Keeping only features associated with first occurences\n");
      }
   }
   else{
      $hExons{$sName}{distinct_exons}{"$sChr|$nStart|$nEnd|$sStrand"}++;
      $hExons{$sName}{first}{start} = $nStart;
      $hExons{$sName}{first}{end}   = $nEnd;
      $hExons{$sName}{first}{score} = $nScore;
      $hExons{$sName}{first}{frame} = $nFrame;
      $hExons{$sName}{last}{start}  = $nStart;
      $hExons{$sName}{last}{end}    = $nEnd;
      $hExons{$sName}{last}{score}  = $nScore;
      $hExons{$sName}{last}{frame}  = $nFrame;
      $hExons{$sName}{chr}     = $sChr;
      $hExons{$sName}{src}     = $sSrc;
      $hExons{$sName}{feature} = $sFeature;
      $hExons{$sName}{strand}  = $sStrand;
      push @aOrder, $sName;
   }
}
close IN;

foreach my $sName (@aOrder){
   unless($flSkipMono and ( scalar(keys %{$hExons{$sName}{distinct_exons}}) == 1 ) ){
      print join("\t", $hExons{$sName}{chr}, $hExons{$sName}{src}, $hExons{$sName}{feature}, $hExons{$sName}{first}{start}, $hExons{$sName}{first}{end},  $hExons{$sName}{first}{score}, $hExons{$sName}{strand}, $hExons{$sName}{first}{frame}, "${sName}|first" ), "\n";
      print join("\t", $hExons{$sName}{chr}, $hExons{$sName}{src}, $hExons{$sName}{feature}, $hExons{$sName}{last}{start}, $hExons{$sName}{last}{end},  $hExons{$sName}{last}{score}, $hExons{$sName}{strand}, $hExons{$sName}{last}{frame}, "${sName}|last" ), "\n";
   }
}
