#!/usr/bin/env perl

# Convert a .wig file into an Affymetrix binary .bar file

# GLOBALS
$ENV{JAVA}    ||= 'java';                                              # Path to java binary
$ENV{WIG2BAR} ||= '/home/harm/hugheslab/opt/USeq_4.1/Apps/Wig2Bar';    # Useq Wig2Bar jar file

# MODULES
use warnings;
use strict;
use Getopt::Long;

# GET ARGUMENTS
my $sHelp       = 0;
my $sInput      = 0;
my $sSeqVersion = '';
my $sSkip       = 0;
GetOptions("help!"       => \$sHelp,
           "file:s"      => \$sInput,
           "version:s"   => \$sSeqVersion,
           "skip:s"      => \$sSkip);


# PRINT HELP
$sHelp = 1 unless($sInput and $sSeqVersion);
if ($sHelp) {
    my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
    die <<HELP

    $sScriptName -f <wigfile> -v <genome-version> -s <skip>
    
    Convert variable step and fixed step wig files into
    Affymetrix bar files, one per chromosome

    Arguments:
    -f <string>
      The full path directory/file name for your xxx.wig(Var) file(s).
    -v <string>
      Genome version (ie H_sapiens_Mar_2006), get from UCSC Browser,
      http://genome.ucsc.edu/FAQ/FAQreleases
    -s <string>
      Skip wig lines with designated value/score.
    -help
      This help message
      
HELP
}


###########
## START ##
###########

# Do the wrap around
my $sArguments = "-f $sInput -s $sSkip ";
$sArguments   .= "-v $sSeqVersion " if ($sSeqVersion);
system("$ENV{JAVA} -jar $ENV{WIG2BAR} $sArguments");
