#!/usr/bin/env perl

# Modified from a bowtie script written by Ben Langmead

# MODULES
use strict;
use warnings;
use Getopt::Long;
use IO::File;
use IO::Zlib qw(:gzip_external 1);

# GET PARAMETERS
my $sHelp    = 0;
my $unpaired = "";
my $mate1    = "";
my $mate2    = "";
my $shuffle  = 0;
GetOptions ("u=s"      => \$unpaired,
            "1=s"      => \$mate1,
            "2=s"      => \$mate2,
            "shuffle!" => \$shuffle,
            "help!"    => \$sHelp) || die "Bad option";

# PRINT HELP
$sHelp = 1 unless($unpaired or ($mate1 and $mate2));
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
    
   Convert a pair of FASTQ files to a single file in a one-pair-per-line
   format, where each line has these five tab-delimited fields:

   1. Pair name
   2. Sequence of mate 1
   3. Qualities of mate 1
   4. Sequence of mate 2
   5. Qualities of mate 2

   This script comes in handy if (a) you'd just like to store your
   paired-end data in a less awkward format than the usual pair of
   parallel FASTQ files, or (b) you'd like to use Bowtie with gzipped
   paired-end input without unzipping it first.  In that case, you can
   pipe the output of this script (which handles gzipped inputs by
   piping them through gzip -dc) to Bowtie and use '--12 -'.
	
   Note that this script can also handle unpaired input (with -u), which
   Bowtie handles approrpaitely in --12 mode even when it's intermingled
   with paired-end input.
   
   Arguments:
    -u <string>
      Name of file(s) with unpaired reads. Multiple file names
      should be comma-separated
    -1 <string>
      Name of file(s) containing the first set of paired reads.
    -2 <string>
      Name of file(s) containing the matched paired reads.
    -shuffle
      Shuffle the read output
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

my @output;
my @unpaireds = split(/,/, $unpaired);

for my $f (@unpaireds) {
   my ($fhUnp);
   if($f =~ /\.gz$/) {
      $fhUnp = new IO::Zlib($f, "rb") or die "Error: can't open '$f': $!\n";
   }
   else{
      $fhUnp = new IO::File($f, "r") or die "Error: can't open '$f': $!\n";
   }
   while(<$fhUnp>) {
      chomp;
      my $name = $_;
      $name = substr($name, 1);
      my $seq = <$fhUnp>;
      chomp($seq);
      my $name2 = <$fhUnp>;
      my $qual = <$fhUnp>;
      chomp($qual);
      if ($shuffle){
	 push @output, "$name\t$seq\t$qual\n";
      }
      else{
	 print "$name\t$seq\t$qual\n";
      }
   }
   close($fhUnp);
}

my @mate1s = split(/,/, $mate1);
my @mate2s = split(/,/, $mate2);

for(my $i = 0; $i <= $#mate1s; $i++) {
   my ($fhM1, $fhM2);
   if($mate1s[$i] =~ /\.gz$/) {
      $fhM1 = new IO::Zlib($mate1s[$i], "rb") or die "Error: can't open '$mate1s[$i]': $!\n";
   }
   else{
      $fhM1 = new IO::File($mate1s[$i], "r") or die "Error: can't open '$mate1s[$i]': $!\n";
   }
   if($mate2s[$i] =~ /\.gz$/) {
      $fhM2 = new IO::Zlib($mate2s[$i], "rb") or die "Error: can't open '$mate2s[$i]': $!\n";
   }
   else{
      $fhM2 = new IO::File($mate2s[$i], "r") or die "Error: can't open '$mate2s[$i]': $!\n";
   }
   
   while(<$fhM1>) {
      my $name1 = $_;
      chomp($name1);
      $name1 = ($name1 =~ / /) ? substr((split(/ /,$name1))[0],1) : substr($name1, 1, -2);
      my $name2 = <$fhM2>;
      chomp($name2);
      $name2 = ($name2 =~ / /) ? substr((split(/ /,$name2))[0],1) : substr($name2, 1, -2);
      my $seq1 = <$fhM1>;
      chomp($seq1);
      my $seq2 = <$fhM2>;
      chomp($seq2);
      my $tmp = <$fhM1>;
      $tmp = <$fhM2>;
      my $qual1 = <$fhM1>;
      chomp($qual1);
      my $qual2 = <$fhM2>;
      chomp($qual2);
      die "Error: mismatch between paired reads ($name1 vs $name2)\n" unless($name1 eq $name2);
      if ($shuffle){
	 push @output, "$name1\t$seq1\t$qual1\t$seq2\t$qual2\n";
      }
      else{
	 print "$name1\t$seq1\t$qual1\t$seq2\t$qual2\n";
      }
   }
   close($fhM1);
   close($fhM2);
}

if($shuffle) {
   @output = shuffle(@output);
   print join("", @output);
}


#################
## SUBROUTINES ##
#################

# Courtesy: http://www.perlmonks.org/?node_id=625977
sub shuffle {
   my @a = \(@_);
   my $n;
   my $i = @_;
   map {
      $n = rand($i--);
      (${$a[$n]}, $a[$n] = $a[$i])[0];
   } @_;
}
