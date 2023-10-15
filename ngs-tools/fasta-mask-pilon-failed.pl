#!/usr/bin/perl

# 08.03.2020 11:20:39 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp         = 0;
my $nBasesPerLine = 100;
my $sVcfFile      = "";
my $sFastaFile    = "";
GetOptions("help!"   => \$sHelp,
           "bases:i" => \$nBasesPerLine,
           "vcf:s"   => \$sVcfFile,
           "fasta:s" => \$sFastaFile);

# PRINT HELP
$sHelp = 1 unless($sVcfFile and $sFastaFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
   
   Arguments:
    -f --fasta <string>
      Pilon-corrected fasta input file
    -v --vcf <string>
      Pilon VCF input file
    -b --bases <integer>
      Number of bases per line in output file
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# READ VCF PASS fields
my %hVCF;
open VCF, "<$sVcfFile" or die "Error: can't read the VCF file: $!\n";
while (<VCF>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sChr, $nPos, $sID, $sRef, $sAlt, $sQual, $sFilter, $sInfo, $sFormat, $sSample) = split /\t/;
   $nPos--;
   $hVCF{$sChr}{$nPos}++ if ($sFilter eq "PASS");
}
close VCF;


# Process the FASTA
my $sFastaHeader = '';
my $sFastaSeq    = '';
my $flRevSeq     = 0;
open INPUT, "<$sFastaFile" or die "Error: can't read the fasta file\n";
while (<INPUT>){
   if (/^>/ or eof){
      if (eof){
         die "Error: file ends in fasta header without sequence\n" if (/^>/);
         $sFastaSeq .= $_;
      }
      if ($sFastaHeader){
         $sFastaHeader =~ s/^>//;
         $sFastaHeader =~ s/ .*$//;
         $sFastaHeader =~ s/_pilon$//;
         $sFastaHeader =~ s/[\n\r]+//g;
         $sFastaSeq  =~ s/\s//g;
         $sFastaSeq  =~ s/[\n\r]+//g;
         
         # Mark non-pass sequences as 'N'
         my $sOutSeq = "";
         foreach (my $i=0; $i < length($sFastaSeq) ; $i++){
            $sOutSeq .= exists( $hVCF{$sFastaHeader}{$i}) ? substr( $sFastaSeq, $i, 1) : "N";
         }
         
         $sOutSeq =~ s/.{$nBasesPerLine}/$&\n/sg;
         $sOutSeq =~ s/\n+$//;
         $sOutSeq .= "\n";
         print join("", ">$sFastaHeader\n", $sOutSeq);
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
