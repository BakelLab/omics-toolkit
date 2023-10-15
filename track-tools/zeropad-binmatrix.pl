#!/usr/bin/env perl

# 02.08.2013 20:06:54 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
GetOptions("help!"   => \$sHelp);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName <matrix> <genelist>
    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

my %hLUT;
my $nF = 0;
open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]'";
while (<IN>){
   print;
   my (@line) = split(/\t/, $_, -1);
   $nF = scalar(@line)>$nF ? scalar(@line) : $nF;
   $hLUT{$line[0]}++;
}
close IN;

print $nF, "\n";

open LIST, $ARGV[1] or die "Error: can't open '$ARGV[1]'";
while (<LIST>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   unless(exists $hLUT{$_}){
      print $_, "\t0"x$nF, "\n";
   }
}
close LIST;
