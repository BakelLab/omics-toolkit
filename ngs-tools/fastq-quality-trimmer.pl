#!/usr/bin/env perl

# 02.09.2010 16:58:25 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use IO::File;
use IO::Zlib qw(:gzip_external 1);

# GET PARAMETERS
my $sHelp       = 0;
my $nMinLen     = 0;
my $nMaxLen     = 0;
my $nQualThresh = 28;
my $nMaxFail    = 10;
my $nMinQual    = 12;
my $sOutput     = '';
my $sEncoding   = 'auto';
my $flNoGap     = 0;
my $flDrop      = 0;
GetOptions("help!"           => \$sHelp,
           "encoding:s"      => \$sEncoding,
           "min_len:i"       => \$nMinLen,
           "max_len:i"       => \$nMaxLen,
           "max_fail_qual:n" => \$nQualThresh,
           "max_fail:i"      => \$nMaxFail,
           "min_qual:i"      => \$nMinQual,
           "output:s"        => \$sOutput,
           "nogap!"          => \$flNoGap,
           "drop!"           => \$flDrop);
my $sIn = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sIn);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [options] <fastq-file>
   
   Process single-end fastq files to remove low-quality read ends. 
   This script assumes an illumina-type chemistry where the phred quality
   of bases progressively declines in a 5' to 3' direction. Output files
   will always be gzipped to save space.
    
   Options:
    -min_qual <integer>
      Minimum phred base quality threshold. Ends of reads will always be cut
      at this base quality. 
      Default: $nMinQual
    -max_fail_qual <integer>
      Reads will be trimmed at the end after encountering 'max_fail' bases 
      below this threshold, counting from start. Trailing bases below the quality
      threshold are also removed.
      Default: $nQualThresh
    -max_fail <integer>
      The number of bases allowed to fall below the 'max_fail_qual' threshold.
      Default: $nMaxFail
    -min_len <integer>
      Minimum required read length after trimming. Trimmed reads shorter than
      this length will not be reported. 
      Default: no limit
    -max_len <integer>
      Maximum read length. All reads will be trimmed to this length at the end.
      Default: no limit
    -output <string>
      Optional output file prefix  
    -nogap
      Remove reads with gaps ('.' or 'N')
    -drop
      Drop rather than filter reads that don't meet filter criteria
    -encoding
      Fastq file encoding (phred33, phred64 or auto). Default: auto
      
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check arguments
die "Error: max_fail_qual must be an integer\n"     unless ($nQualThresh =~ /^\d+$/);
die "Error: max_fail must be between 0 and 1\n"   unless ($nMaxFail =~ /^\d+$/);
die "Error: min_len must be an integer\n"         unless ($nMinLen =~ /^\d+$/);
die "Error: max_len must be an integer\n"         unless ($nMaxLen =~ /^\d+$/);
die "Error: max_len must be greater than min_len" if($nMaxLen and ($nMaxLen<$nMinLen));

# Set file encoding
my $nEncoding = 0;
if (lc($sEncoding) eq "auto"){
   $nEncoding = auto_fastq_encoding($sIn);
   die "Error: could not determine fastq encoding automatically, please specify manually\n" unless($nEncoding);
   warn "Encoding set to phred${nEncoding}\n";
}
elsif (lc($sEncoding) eq "phred33"){
   $nEncoding = 33;
}
elsif (lc($sEncoding) eq "phred64"){
   $nEncoding = 64;
}
else{
   die "Error: unknown fastq encoding specified, use one of 'phred33', 'phred64' or 'auto'.\n";
}

# Open files
my @asSuffixes = ('.fq','.fq.gz','.FQ','.FQ.GZ','.fastq','.fastq.gz','.FASTQ','.FASTQ.GZ');
my $sOut  = $sOutput ? "${sOutput}.fastq.gz"   : join('', (fileparse($sIn, @asSuffixes))[0], '.trimmed.fastq.gz');
my $sStat = $sOutput ? "${sOutput}.fastq.stat" : join('', (fileparse($sIn, @asSuffixes))[0], '.trimmed.fastq.stat');
my $fhIn  = gzopen($sIn);
my $fhOut = new IO::Zlib($sOut, "wb6") or die "Error: can't open '$sOut': $!\n";

# Process the reads
my %hStats;
my $nCount = 0;
while (<$fhIn>){
   my $sName = $_;
   my $sSeq  = <$fhIn>;
   my $sTmp  = <$fhIn>;
   my $sQual = <$fhIn>;
   $sName =~ s/[\n\r]+$//;
   $sSeq  =~ s/[\n\r]+$//;
   $sQual =~ s/[\n\r]+$//;
   
   # Set start length
   my $nStartLength = length($sSeq);
   
   # Only include reads without gaps
   next if ($sSeq =~ /\.|N/ and $flNoGap);
   
   # Trim size to maxlength
   die "Error: length mismatch between sequence and quality strings\n" unless(length($sSeq) == length($sQual));
   if ($nMaxLen){
      $sSeq = substr($sSeq,   0, $nMaxLen) if (length($sSeq)  > $nMaxLen);
      $sQual = substr($sQual, 0, $nMaxLen) if (length($sQual) > $nMaxLen);
   }
   
   # Process Read quality
   my $nLength    = 0;
   if ($nQualThresh){
      my @asSegQual = split //, $sQual;
      my $nConsCount = 0;
      my $nFailCount = 0;
      for(my $i=0 ; $i<@asSegQual ; $i++){
         $nLength  = $i+1;
         my $nSegQual = ord($asSegQual[$i]) - $nEncoding;
         if ($nSegQual < $nQualThresh){
            $nConsCount++;
            if ($nSegQual <= $nMinQual){
               last;
            }
            else{
               $nFailCount++;
               if ($nFailCount > $nMaxFail){
                  last;
               }
            }
         }
         else{
            $nConsCount = 0
         }
      }
      $nLength = $nLength - $nConsCount;
   }
   else{
      $nLength = length($sSeq);
   }
   next if ($nLength < $nMinLen);
   next if ($nLength != $nStartLength and $flDrop);
   
   # Write the fastq output
   print $fhOut join('', $sName, "\n", substr($sSeq, 0, $nLength), "\n+\n", substr($sQual, 0, $nLength)), "\n";
   
   # Keep stats
   $hStats{length} += $nLength;
   $hStats{distrib}{$nLength}++;
   $nCount++;
}
close $fhIn;
close $fhOut;


# Output the stats
open STAT, ">$sStat" or die "Error: can't open '$sStat': $!\n";
print STAT "#Read count: $nCount\n";
print STAT "#Base count: $hStats{length}\n";
for my $nLength (sort {$a <=> $b} keys(%{$hStats{distrib}})){print STAT "$nLength\t$hStats{distrib}{$nLength}\n"}
close STAT;


###############
# SUBROUTINES #
###############

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


# auto_fastq_encoding
#
# Try to automatically determine the encoding of a fastq file
# Returns 33 or 64 for phred33 or phred64 files respectively, or 0 if automatic scan failed
sub auto_fastq_encoding {
   my ($sFastq) = @_;
   my $nReturn  = 0;
   
   my $fhFastq = gzopen($sFastq);
   while (<$fhFastq>){
      my $nPhred33count = 0;
      my $nPhred64count = 0;
      my $sName = $_;
      my $sSeq  = <$fhFastq>;
      my $sTmp  = <$fhFastq>;
      my $sQual = <$fhFastq>;
      $sQual =~ s/[\n\r]+$//;
      my @asSegQual = split //, $sQual;
      foreach my $sSegQual (@asSegQual){
         my $nSegQual = ord($sSegQual);
         if ($nSegQual <59){
            $nPhred33count++;
         }
         if ($nSegQual >74){
            $nPhred64count++;
            
         }
      }
      if ($nPhred33count > 2 && $nPhred64count == 0){
         $nReturn = 33;
         last;
      }
      elsif ($nPhred33count == 0 && $nPhred64count > 2){
         $nReturn = 64;
         last;
      }
      else{
         last if $. > 2000;
      }
   }
   close $fhFastq;
   return $nReturn;
}

