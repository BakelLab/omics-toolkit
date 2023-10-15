#!/usr/bin/env perl

# 09.09.2011 11:15:48 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
my $sHeader      = '';
GetOptions("help!"    => \$sHelp,
           "header:s" => \$sHeader);
my $sInput = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sInput and $sHeader);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -h header <vcf-file>
    
    -header <string>
      Comma-separated header info for the multi-sample calls
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########



# Parse the header
my @asHeader = split /,/, $sHeader;

# Now process the vcf file
my %hCalls;
open VCF, $sInput or die "Error: can't open '$sInput': $!\n";
while (<VCF>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my ($sChr, $nPos, $sID, $sRef, $sAlt, $sQual, $sFilter, $sInfo, $sFormat, @asGenotypes) = split /\t/, $_, -1;
   die "Error: number of header fields don't match the number of genotypes\n" unless(@asHeader == @asGenotypes);
   
   # skip multicalls
   next if($sAlt =~ /,/);
   
   # Clean up the genotype calls
   for (my $i=0 ; $i<@asGenotypes ; $i++){
      $asGenotypes[$i] =~ s/\:.*$//;
   }
   
   # Get the heterozygosity calls in each sample
   my @asRefAlleles = ($sRef, $sAlt);
   for (my $i=0 ; $i<@asGenotypes ; $i++){
      my ($a, $b) = split /\/|\|/, $asGenotypes[$i];
      my $sSNVcall = join('', $asRefAlleles[$a], $asRefAlleles[$b]);
      $sSNVcall =~ s/,//g;
      $hCalls{$asHeader[$i]} .= get_IUPAC($sSNVcall);
   }   
}
close VCF;

# Print fasta output
foreach my $sKey (keys %hCalls){
   print ">$sKey\n$hCalls{$sKey}\n";
}


###############
# SUBROUTINES #
###############


# get_IUPAC
#
# Get IUPAC code according to bases in a string
sub get_IUPAC {
   my ($sString) = @_;
   
   #IUPAC CODES
   my %hIUPAC = ( "A"    => "A",
                  "T"    => "T",
                  "C"    => "C",
                  "G"    => "G",
                  "AC"   => "M",
                  "AG"   => "R",
                  "AT"   => "W",
                  "CG"   => "S",
                  "CT"   => "Y",
                  "GT"   => "K",
                  "ACG"  => "V",
                  "ACT"  => "H",
                  "AGT"  => "D",
                  "CGT"  => "B",
                  "ACGT" => "N");
   
   # Parse string
   my %hChars;
   my @asChars = split //, $sString;
   foreach my $sChar (@asChars) {$hChars{$sChar}++;}
   my $sMatch = join('', sort(keys %hChars));
   if (exists $hIUPAC{$sMatch}){
      return $hIUPAC{$sMatch};
   }
   else{
      die "Error: non-base characters found in IUPAC string conversion: '$sMatch'\n";
   }
}


