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
my $sHelp  = 0;
my $flAll  = 0;
GetOptions("help!"    => \$sHelp,
           "all!"     => \$flAll);

# PRINT HELP
$sHelp = 1 unless(@ARGV);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [options] <qseq-file1> .. <qseq-fileN>
   
   Converts illumina qseq files to fastq format. By default it
   only keeps reads that passed Illumina QC filters.
    
   Options:
    -all
      Also include reads that did not pass QC filters
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

foreach my $sIn (@ARGV){
   # Open file depending on extension
   my ($fhIn);
   if($sIn =~ /\.gz$/) {
      $fhIn = new IO::Zlib($sIn, "rb") or die "Error: can't open '$sIn': $!\n";
   }
   else{
      $fhIn = new IO::File($sIn, "r") or die "Error: can't open '$sIn': $!\n";
   }
   
   # Process the reads
   while (<$fhIn>){
      s/[\n\r]+$//g;
      my @asRead   = split /\t/;
      next unless ($asRead[$#asRead] or $flAll);
      my $sRead_id = join(':', @asRead[0,2,3,4,5,7]);
      print join('', '@', $sRead_id, "\n", $asRead[8], "\n+\n", $asRead[9]), "\n";
   }
   close $fhIn;
}
