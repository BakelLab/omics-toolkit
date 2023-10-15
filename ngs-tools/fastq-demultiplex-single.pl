#!/usr/bin/env perl

# 01.12.2011 15:25:09 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use String::Approx 'amatch';
use File::Basename;
use IO::File;
use IO::Zlib qw(:gzip_external 1);

# GET PARAMETERS
my $sHelp        = 0;
my $sFastqFile   = "";
my $sBarcodeFile = "";
my $nMismatches  = 1;
GetOptions("help!"      => \$sHelp,
           "fastq:s"    => \$sFastqFile,
           "barcode:s"  => \$sBarcodeFile,
           "mismatch:i" => \$nMismatches);

# PRINT HELP
$sHelp = 1 unless($sFastqFile and $sBarcodeFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
   
   Script to demultiplex a fastq file that has barcodes at the start
   of each read. The barcodes should be provided in a separate file with
   three columns: barcode, sampleID, left_trim. The last column
   specifies how much of the beginning of the read should be trimmed off.
   Demultiplexed files are gzipped by default to conserve space.
   
   Arguments:
    -f --fastq <string>
      Fastq file to demultiplex
    -b --barcode <string>
      Tab-delimited file with barcodes (barcode, sampleID, left_trim)
    -m --mismatch <integer>
      Allowed number of mismatches. Default: $nMismatches
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Format the output file prefix
my ($sOutputPrefix,$sPath,$sSuffix) = fileparse($sFastqFile,qr{\.fastq\..*$});

# Read the barcode file and initialize output filehandles
my ($rBarcodes, $nBarcodeLength)  = init_barcodes($sBarcodeFile, $sOutputPrefix);
my @asBarcodes = keys(%$rBarcodes);

# Open container for unmatched barcodes
my $fhUnknown = new IO::Zlib("unmatched_${sOutputPrefix}.fastq.gz", "wb6") or die "Error: can't write to '${sOutputPrefix}_unmatched.fastq.gz': $!\n";

# Process the sequences
my $fhIn = gzopen($sFastqFile);
while (<$fhIn>){
   my $sName = $_;
   my $sSeq  = <$fhIn>;
   my $sTmp  = <$fhIn>;
   my $sQual = <$fhIn>;
   $sName =~ s/[\n\r]+$//;
   $sTmp  =~ s/[\n\r]+$//;
   $sSeq  =~ s/[\n\r]+$//;
   $sQual =~ s/[\n\r]+$//;
   die "Error: '$sFastqFile' is not a fastq file\n" unless( ($sName =~ /^@/) and ($sTmp =~ /^\+/) );
   if (length($sSeq) > $nBarcodeLength){
      my @asMatches = amatch(substr($sSeq,0,$nBarcodeLength), ["S$nMismatches"], @asBarcodes);
      if ( (@asMatches == 1) and (length($sSeq) > $rBarcodes->{$asMatches[0]}{trimleft}) ){
         $sSeq  = substr($sSeq,$rBarcodes->{$asMatches[0]}{trimleft});
         $sQual = substr($sQual,$rBarcodes->{$asMatches[0]}{trimleft});
         $rBarcodes->{$asMatches[0]}{filehandle}->print("$sName\n$sSeq\n$sTmp\n$sQual\n");
      }
      else{
         $fhUnknown->print("$sName\n$sSeq\n$sTmp\n$sQual\n");
      }
   }
   else{
      $fhUnknown->print("$sName\n$sSeq\n$sTmp\n$sQual\n");
   }
}

# Close the output filehandles
foreach my $sBarcode (@asBarcodes){
   $rBarcodes->{$sBarcode}{filehandle}->close();
}
close $fhUnknown;



#################
## SUBROUTINES ##
#################

# gzopen
#
# Wrapper to open a file in uncompressed or compressed format
sub gzopen {
   my ($sFile) = @_;
   my $fhFile;
   if ($sFile =~ /\.gz$/){
      $fhFile = new IO::Zlib($sFile, "rb") or die "Error: can't open '$sFile': $!\n";
   }
   else{
      $fhFile = new IO::File($sFile, "r") or die "Error: can't open '$sFile': $!\n";
   }
   return $fhFile;
}


# init_barcodes($barcodefile)
#
# Reads the barcode file into a hash
sub init_barcodes {
   my ($sBarcodeFile, $sOutputPrefix) = @_;
   my %hBarcodes;
   my $nMaxLen = 0;
   open BARCODE, $sBarcodeFile or die "Error: can't open '$sBarcodeFile': $!\n";
   while (<BARCODE>){
      next if (/^\s*$/);
      next if (/^ *#/);
      s/[\n\r]+$//;
      my($sBarcode, $sSampleID, $nTrimLeft) = split /\t/;
      $sSampleID =~ s/^ +//;
      $sSampleID =~ s/ +$//;
      $sSampleID =~ s/ /_/g;
      my $sOutFile  = join("", $sSampleID, "_", $sOutputPrefix, ".fastq.gz");
      die "Error: non-base characters found in barcode '$sBarcode'\n" unless($sBarcode =~ /^[ACGTacgt]+$/);
      die "Error: '$nTrimLeft' is not an integer in '$sBarcodeFile', line $.\n" unless($nTrimLeft =~ /^\d+$/);
      die "Error: duplicate barcode '$sBarcode' found in barcode file\n" if (exists $hBarcodes{$sBarcode});
      $hBarcodes{$sBarcode}{trimleft}   = $nTrimLeft;
      $hBarcodes{$sBarcode}{outfile}    = $sOutFile;
      $nMaxLen = length($sBarcode) if (length($sBarcode) > $nMaxLen);
   }
   close BARCODE;
   
   # Now create filehandles for all the barcodes and return the barcode hash
   foreach my $sBarcode (keys %hBarcodes){
      $hBarcodes{$sBarcode}{filehandle} = new IO::Zlib($hBarcodes{$sBarcode}{outfile}, "wb6") or die "Error: can't write to '$hBarcodes{$sBarcode}{outfile}': $!\n";
   }
   return(\%hBarcodes, $nMaxLen);
}

