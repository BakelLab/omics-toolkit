#!/usr/bin/env perl

# 07.09.2010 20:00:48 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp         = 0;
my $nMinLength    = 0;
my $nMaxLength    = 0;
GetOptions("help!"       => \$sHelp,
           "minlength:i" => \$nMinLength,
           "maxlength:i" => \$nMaxLength);

# PRINT HELP
$sHelp = 1 unless(@ARGV);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [options] <fasta-file1> .. <fasta-fileN>
   
   Extract sequences from a multi-fasta file that meet
   certain length thresholds.
    
   Options:
    -minlength <integer>
      Minumum length of fasta sequence
      Default: 0 (no filtering)
    -maxlength <integer>
      Maximum length of fasta sequence
      Default: 0 (no filtering)
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

foreach my $sInput (@ARGV){
   my $sFastaHeader = '';
   my $sFastaSeq    = '';
   my $flRevSeq     = 0;
   open INPUT, "<$sInput" or die "Error: can't read the fasta file\n";
   while (<INPUT>){
      if (/^>/ or eof){
         if (eof){
            die "Error: file ends in fasta header without sequence\n" if (/^>/);
            $sFastaSeq .= $_;
         }
         if ($sFastaHeader){
            my $flPrint = 1;
            my $sFastaLengthCheck = $sFastaSeq;
            $sFastaLengthCheck  =~ s/\s//g;
            $sFastaLengthCheck  =~ s/[\n\r]+//g;
            $flPrint = 0 if ($nMinLength and (length($sFastaLengthCheck) < $nMinLength));
            $flPrint = 0 if ($nMaxLength and (length($sFastaLengthCheck) > $nMaxLength));
            if($flPrint){
               print join("", $sFastaHeader, $sFastaSeq);
            }
         }
         $sFastaHeader = $_;
         $sFastaSeq    = "";
      }
      else{
         next if (/^\s*$/);
         next if (/^ *#/);
         $sFastaSeq .= $_ if ($sFastaHeader);
      }
   }
   close INPUT;
}
