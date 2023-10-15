#!/usr/bin/env perl

# 01.11.2012 12:21:24 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use File::Basename;

# TEMP DIR LOCATION
$ENV{TMPSHMDIR} ||= '/dev/shm';
$ENV{PHYMMBL}   ||= '/projects/vanbah01a/opt/phymmbl';

# GET PARAMETERS
my $sHelp        = 0;
GetOptions("help!"   => \$sHelp);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName <file.fa>
   
   Wraps the phymmbl scoreReads script to run on /dev/shm
   to speed things up and to avoid temp files from mixing.
    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Get current dir
my $sStartDir = getcwd();
my $sBasename = basename($ARGV[0]);

# Set up all the necessary files in the temp dir
system("cp $ENV{PHYMMBL}/scoreReads.pl $ENV{TMPSHMDIR}")  == 0 or die "Error: could not copy scoring script: $!\n";
system("cp $ARGV[0] $ENV{TMPSHMDIR}")   == 0 or die "Error: could not copy fasta file: $!\n";
system("cp -R $ENV{PHYMMBL}/.blastData $ENV{TMPSHMDIR}") == 0 or die "Error: could not copy blast data: $!\n";
system("ln -s $ENV{PHYMMBL}/.genomeData $ENV{TMPSHMDIR}/.genomeData") == 0 or die "Error: could not symlink genome data: $!\n";
system("ln -s $ENV{PHYMMBL}/.scripts $ENV{TMPSHMDIR}/.scripts") == 0 or die "Error: could not symlink script data: $!\n";
system("ln -s $ENV{PHYMMBL}/.taxonomyData $ENV{TMPSHMDIR}/.taxonomyData") == 0 or die "Error: could not symlink taxonomy data: $!\n";
system("mkdir -p $ENV{TMPSHMDIR}/.logs") == 0 or die "Error: could not create logs dir: %!\n";

# Run the stuff in the temp dir
chdir($ENV{TMPSHMDIR});
system("./scoreReads.pl $sBasename")  == 0 or die "Error: could not run scoreReads.pl: %!\n";

# Copy the result files
system("cp results.* $sStartDir") == 0 or die "Error: could not copy final results\n";

# Clean stuff up
system("rm -f scoreReads.pl $sBasename .genomeData .scripts .taxonomyData");
system("rm -rf .blastData .logs");


