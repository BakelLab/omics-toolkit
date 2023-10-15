#!/usr/bin/env perl

# Convert a .sgr file into an Affymetrix binary .bar file

# MODULES
use strict;
use Getopt::Long;
use Convert::Binary::C;
use File::Basename;

# GET ARGUMENTS
my $sHelp       = 0;
my $sSeqGroup   = 0;
my $sSeqVersion = '';
my $sBarFile    = '';
GetOptions("help!"       => \$sHelp,
           "group:s"     => \$sSeqGroup,
           "version:s"   => \$sSeqVersion,
           "bar:s"       => \$sBarFile);


# PRINT HELP
$sHelp = 1 unless(@ARGV);
if ($sHelp) {
    die <<HELP

    sgr2bar.pl [-g group -v version] <sgr-file>

    Arguments:
    -group <string>
      The sequence group name, usually an abbreviation of the organism name
      examples: 'Mm' or 'Hs' for human and mouse respectively
    -version <string>
      The sequence version, i.e. the genome sequence version
      examples: NCBI36
    -bar <string>
      The name of the output file. If this argument is omitted, the name
      will be based on the input file.
    -help
      This help message
      
HELP
}


###########
## START ##
###########

# First we need to get a listing of all chromosomes and the number of probes per chromosome
my $sSgrFile         = shift @ARGV;
my $rhChrProbeCounts = &get_probe_counts($sSgrFile);

# Start printing the binary output
my ($sName, $sPath, $sSuffix) = fileparse($sSgrFile, qr/\.[^.]*/);
$sName = $sBarFile ? $sBarFile : "$sName.bar";
open BIN, ">$sName" and binmode(BIN);

# First print the magic number to identify the filetype
print BIN pack("a4", 'barr');
print BIN pack("a", "\r");
print BIN pack("a","\n");
print BIN pack("a", "\032");
print BIN pack("a", "\n");

# The bar file header
my $c = Convert::Binary::C->new(ByteOrder => 'BigEndian');
print BIN $c->pack('float', 2.0);                                 # Bar file version
print BIN $c->pack('int', scalar(keys(%$rhChrProbeCounts)));      # No. of sequences (chromosomes) in the file
print BIN $c->pack('int', 2);                                     # Number of data columns
print BIN $c->pack('int', 2);                                     # Field type 1 'integer' (position)
print BIN $c->pack('int', 1);                                     # Field type 2 'float'   (value)
print BIN $c->pack('int', 0);                                     # The number of parameter-value pairs

# Now start processing the measurements in the file again, this time packing into binary output
# Note that the input file is sorted by chromosome and position before processing
my $sCurrentChr = '';
open IN, "sort -t '\t' -k1,1b -k2,2n $sSgrFile |" or die "Can't open sgr file: $!\n";
while(<IN>){
   next if /^\s*$/;
   next if /^\s*#/;
   s/[\n\r]//g;
   s/ //g;
   my ($sChr, $nPos, $nVal) = split /\t/;
   if ($sChr ne $sCurrentChr){                                    # New chromosome; print sequence header
      my $nSeqNameLength = length($sChr);
      print BIN $c->pack('int', $nSeqNameLength);                 # Sequence name length
      print BIN pack("a$nSeqNameLength", $sChr);                  # Sequence name
      
      if ($sSeqGroup){
         $sSeqGroup =~ s/ /_/g;
         my $nSeqGroupLength = length($sSeqGroup);
         print BIN $c->pack('int', $nSeqGroupLength);                # Sequence group name length
         print BIN pack("a$nSeqGroupLength", $sSeqGroup);            # Sequence group name
      }
      else{
         print BIN $c->pack('int', 0);
      }
      
      if ($sSeqVersion){
         $sSeqVersion =~ s/ /_/g;
         my $nSeqVersionLength = length($sSeqVersion);
         print BIN $c->pack('int', $nSeqVersionLength);              # Version length
         print BIN pack("a$nSeqVersionLength", $sSeqVersion);        # Version name
      }
      else{
         print BIN $c->pack('int', 0);
      }
      
      print BIN $c->pack('int', 0);                               # no. of par/val pairs
      print BIN $c->pack('int', $rhChrProbeCounts->{$sChr});      # no. of data points
      
      $sCurrentChr = $sChr;
   }
   
   print BIN $c->pack('int', $nPos);
   print BIN $c->pack('float', $nVal);
}

# Close the file
print BIN pack("a4", "END\n");
close BIN;


#################
## SUBROUTINES ##
#################

# get_probe_counts(filename)
#
# Reads the sgr file and returns the number of probes per chromosome in a hash
# Also does some checking along the way
sub get_probe_counts{
   my $sInput = shift @_;
   my %hReturn;
   my $nCount = 1;
   open SGR, $sInput or die "Can't open sgr file: $!\n";
   while(<SGR>){
      next if /^\s*$/;
      next if /^\s*#/;
      s/ //g;
      chomp;
      my @asLine = split /\t/;
      die "Error: line $nCount does not have the required number of fields\n" unless(scalar(@asLine)==3);
      die "Error: no chromosome specified on line $nCount\n" unless ($asLine[0]);
      die "Error: position should be an integer on line $nCount\n" unless ($asLine[1] =~ /^\d+$/);
      die "Error: value should be a float on line $nCount\n" unless ($asLine[2] =~ /^-?\d+\.?\d*e?[+-]?\d*$/);
      $hReturn{$asLine[0]}++;
      $nCount++;
   }
   return \%hReturn;
   close SGR;
}
