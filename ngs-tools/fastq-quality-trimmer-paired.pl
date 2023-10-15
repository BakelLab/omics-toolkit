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
my $sM1         = '';
my $sM2         = '';
my $nMinLen     = 55;
my $nMaxLen     = 0;
my $nQualThresh = 28;
my $nMinQual    = 12;
my $nMaxFail    = 10;
my $sOutput     = '';
my $sEncoding   = 'auto';
my $flNoGap     = 0;
GetOptions("help!"           => \$sHelp,
           "1:s"             => \$sM1,
           "2:s"             => \$sM2,
           "min_len:i"       => \$nMinLen,
           "max_len:i"       => \$nMaxLen,
           "max_fail_qual:n" => \$nQualThresh,
           "max_fail:i"      => \$nMaxFail,
           "min_qual:i"      => \$nMinQual,
           "output:s"        => \$sOutput,
           "encoding:s"      => \$sEncoding,
           "nogap!"          => \$flNoGap);

# PRINT HELP
$sHelp = 1 unless($sM1 and $sM2);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [options] -1 <fastq-file1> -2 <fastq-file2>
    
   Process paired-end fastq files to remove low-quality read ends. 
   This script assumes an illumina-type chemistry where the phred quality
   of bases progressively declines in a 5' to 3' direction. Output files
   will always be gzipped to save space.
   
    -1 <string>
      Read file 1 (fastq or fastq.gz format)
    -2 <string>
      Read file 2 (fastq or fastq.gz format)
   
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
      Drop reads with gaps ('.' or 'N')
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
die "Error: max_fail must be between 0 and 1\n"     unless ($nMaxFail =~ /^\d+$/);
die "Error: min_len must be an integer\n"           unless ($nMinLen =~ /^\d+$/);
die "Error: max_len must be an integer\n"           unless ($nMaxLen =~ /^\d+$/);
die "Error: max_len must be greater than min_len\n" if($nMaxLen and ($nMaxLen<$nMinLen));
die "Error: fastq files are identical\n"            if ($sM1 eq $sM2);

# Set file encoding
my $nEncoding = 0;
if (lc($sEncoding) eq "auto"){
   my $nEncoding1 = auto_fastq_encoding($sM1);
   my $nEncoding2 = auto_fastq_encoding($sM2);
   if ($nEncoding1 and !$nEncoding2){
      $nEncoding = $nEncoding1;
   }
   elsif (!$nEncoding1 and $nEncoding2){
      $nEncoding = $nEncoding2;
   }
   elsif ($nEncoding1 and $nEncoding2){
      die "Error: conflicting fastq encoding ('$nEncoding1' vs '$nEncoding2') for paired files\n" unless ($nEncoding1 == $nEncoding2);
      $nEncoding = $nEncoding1;
   }
   else{
      die "Error: could not determine fastq encoding automatically, please specify manually\n";
   }
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
my $sO1    = $sOutput ? "${sOutput}_1.fastq.gz"   : join('', (fileparse($sM1, @asSuffixes))[0], '.trimmed.fastq.gz');
my $sO2    = $sOutput ? "${sOutput}_2.fastq.gz"   : join('', (fileparse($sM2, @asSuffixes))[0], '.trimmed.fastq.gz');
my $sStat1 = $sOutput ? "${sOutput}_1.fastq.stat" : join('', (fileparse($sM1, @asSuffixes))[0], '.trimmed.fastq.stat');
my $sStat2 = $sOutput ? "${sOutput}_2.fastq.stat" : join('', (fileparse($sM2, @asSuffixes))[0], '.trimmed.fastq.stat');
my $fhM1 = gzopen($sM1);
my $fhM2 = gzopen($sM2);
my $fhO1 = new IO::Zlib($sO1, "wb6") or die "Error: can't open '$sO1': $!\n";
my $fhO2 = new IO::Zlib($sO2, "wb6") or die "Error: can't open '$sO2': $!\n";

# Process the reads
my %hStats1;
my %hStats2;
my $nCount = 0;
while (<$fhM1>){
   my $name1 = $_;
   my $name2 = <$fhM2>;
   my $seq1  = <$fhM1>;
   my $seq2  = <$fhM2>;
   my $tmp1  = <$fhM1>;
   my $tmp2  = <$fhM2>;
   my $qual1 = <$fhM1>;
   my $qual2 = <$fhM2>;
   $name1 =~ s/[\n\r]+$//;
   $name2 =~ s/[\n\r]+$//;
   $name1 =~ s/\s.*$//;
   $name2 =~ s/\s.*$//;
   $seq1  =~ s/[\n\r]+$//;
   $seq2  =~ s/[\n\r]+$//;
   $qual1 =~ s/[\n\r]+$//;
   $qual2 =~ s/[\n\r]+$//;
   $name1 =~ s/_Read1//;
   $name2 =~ s/_Read2//;
   if (length($qual1)>=0 and length($qual2)>=0){
      die "Error: read name '$name1' in file '$sM1' does not conform to the fastq standard\n" unless ($name1 =~ /^\@/);
      die "Error: read name '$name2' in file '$sM2' does not conform to the fastq standard\n" unless ($name2 =~ /^\@/);
      unless(substr($name1,1,-2) eq substr($name2,1,-2)){
         warn "Warning: mismatch between paired reads ($name1 vs $name2), skipping\n" ;
         next;
      }
   }
   else{
      die "Error: Missing paired reads: $name1  <=>  $name2\n";
   }
   
   # Only include reads without gaps
   next if ( (($seq1 =~ /[nN.]/) or ($seq2 =~ /[nN.]/)) and $flNoGap);
   
   # Trim to maximum length
   die "Error: length mismatch between sequence and quality strings for read '$name1'\n" unless(length($seq1) == length($qual1));
   die "Error: length mismatch between sequence and quality strings for read '$name2'\n" unless(length($seq2) == length($qual2));
   if ($nMaxLen){
      $seq1  = substr($seq1, 0, $nMaxLen)  if (length($seq1) > $nMaxLen);
      $qual1 = substr($qual1, 0, $nMaxLen) if (length($qual1) > $nMaxLen);
      $seq2  = substr($seq2, 0, $nMaxLen)  if (length($seq2) > $nMaxLen);
      $qual2 = substr($qual2, 0, $nMaxLen) if (length($qual2) > $nMaxLen);
   }
   
   # Now move on to the more time consuming tests; filter on segment quality
   my $nLength1    = 0;
   my $nLength2    = 0;
   if ($nQualThresh){
      my @asSegQual1 = split //, $qual1;
      my @asSegQual2 = split //, $qual2;

      # Process read1 quality
      my $nConsCount1 = 0;
      my $nFailCount1 = 0;
      for(my $i=0 ; $i<@asSegQual1 ; $i++){
         $nLength1 = $i+1;
         my $nSegQual1 = ord($asSegQual1[$i]) - $nEncoding;
         if ($nSegQual1 < $nQualThresh){
            $nConsCount1++;
            if ($nSegQual1 <= $nMinQual){
               last;
            }
            else{
               $nFailCount1++;
               if ($nFailCount1 > $nMaxFail){
                  last;
               }
            }
         }
         else{
            $nConsCount1 = 0
         }
      }
      $nLength1 = $nLength1 - $nConsCount1;
      next if ($nLength1 < $nMinLen);
      
      # Process read2 quality
      my $nConsCount2 = 0;
      my $nFailCount2 = 0;
      for(my $i=0 ; $i<@asSegQual2 ; $i++){
         $nLength2 = $i+1;
         my $nSegQual2 = ord($asSegQual2[$i]) - $nEncoding;
         if ($nSegQual2 < $nQualThresh){
            $nConsCount2++;
            if ($nSegQual2 <= $nMinQual){
               last;
            }
            else{
               $nFailCount2++;
               if ($nFailCount2 > $nMaxFail){
                  last;
               }
            }
         }
         else{
            $nConsCount2 = 0
         }
      }
      $nLength2 = $nLength2 - $nConsCount2;
      next if ($nLength2 < $nMinLen);
   }
   else{
      $nLength1 = length($seq1);
      $nLength2 = length($seq2);
      next if ($nLength1 < $nMinLen);
      next if ($nLength2 < $nMinLen);
   }
   
   # Write the fastq output
   print $fhO1 join('', $name1, "\n", substr($seq1, 0, $nLength1), "\n+\n", substr($qual1, 0, $nLength1)), "\n";
   print $fhO2 join('', $name2, "\n", substr($seq2, 0, $nLength2), "\n+\n", substr($qual2, 0, $nLength2)), "\n";
   
   # Keep stats
   $hStats1{length} += $nLength1;
   $hStats1{distrib}{$nLength1}++;
   $hStats2{length} += $nLength2;
   $hStats2{distrib}{$nLength2}++;
   $nCount++;
}
close $fhM1;
close $fhM2;
close $fhO1;
close $fhO2;


# Output the stats
open STAT1, ">$sStat1" or die "Error: can't open '$sStat1': $!\n";
print STAT1 "#Read count: $nCount\n";
if (exists $hStats1{length}){
   print STAT1 "#Base count: $hStats1{length}\n";
   for my $nLength (sort {$a <=> $b} keys(%{$hStats1{distrib}})){print STAT1 "$nLength\t$hStats1{distrib}{$nLength}\n"}
}
close STAT1;

open STAT2, ">$sStat2" or die "Error: can't open '$sStat2': $!\n";
print STAT2 "#Read count: $nCount\n";
if (exists $hStats2{length}){
   print STAT2 "#Base count: $hStats2{length}\n";
   for my $nLength (sort {$a <=> $b} keys(%{$hStats2{distrib}})){print STAT2 "$nLength\t$hStats2{distrib}{$nLength}\n"}
}
close STAT2;


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
