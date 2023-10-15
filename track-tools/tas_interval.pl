#!/usr/bin/env perl

# Wrapper around tilecore interval analysis that allows for specification of a
# custom output filename

use strict;
use warnings;
use Getopt::Long;
use Cwd;
use File::Temp qw(tempfile tempdir);
use File::Basename;
use File::Path;
use File::Copy;

# GLOBALS
$ENV{TMPDIR}    ||= "/tmp";                                         # location for tmp file storage
$ENV{TILECORE}  ||= "/home/hugheslab1/projects/bin/tilecore.exe";   # location of tilecore binary
$SIG{'INT'}=$SIG{'HUP'}=$SIG{'ABRT'}=$SIG{'QUIT'}=$SIG{'TRAP'}=$SIG{'STOP'}=\&INTERRUPT;

# GET ARGUMENTS
my $sHelp       = 0;
my $sBarFile    = '';
my $sOutFile    = '';
my $nMaxGap     = 100;
my $nMinRun     = 60;
my $nThreshold  = 0;
my $nLess       = 0;
GetOptions("help!"   => \$sHelp,
           "bar:s"   => \$sBarFile,
           "out:s"   => \$sOutFile,
           "gap:i"   => \$nMaxGap,
           "run:i"   => \$nMinRun,
           "thr:f"   => \$nThreshold,
           "less!"   => \$nLess);


# PRINT HELP
$sHelp = 1 unless($sBarFile);
if ($sHelp) {
    my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
    die <<HELP

    $sScriptName

    Arguments:
    -bar <string>
      The bar file
    -gap <integer>
      max_gap, default: $nMaxGap
    -run <integer> 
      min_run, default: $nMinRun
    -thr <number>
      threshold, default: $nThreshold
    -out
      The name of the output file (optional)
    -less
      If specified, probes below the threshold are selected. The default is to
      select probes scoring above the threshold.
    -help
      This help message
      
HELP
}

###########
## START ##
###########

# Argument checking
$nLess = 1 if ($nLess);

# Get current working dir and create a temporary folder
my $sCwd    = getcwd();
my $sTmpDir = tempdir('tilecore-XXXXX', DIR=>$ENV{TMPDIR}, CLEANUP=>1);
die "Error: Failed to create temporary dir\n" unless( -e $sTmpDir);

# Create a symlink to the bar file in the temporary dir
my $sBarName = basename($sBarFile);
my $sBarLink = $sBarFile =~ /^\// ? $sBarFile : "$sCwd/$sBarFile";
symlink($sBarLink, "$sTmpDir/$sBarName") or die "Error: Could not create symlink to '$sBarFile' file in the temporary folder: $!\n";
die "Error: Could not create symlink to '$sBarFile' file in the temporary folder\n" unless ( -e "$sTmpDir/$sBarName" );

# Run tilecore in the temp dir and when it's done rename the output file and move it
# to the folder from which the script was run
chdir $sTmpDir;
my $sCommand = "/usr/bin/wine $ENV{TILECORE} -interval -bar $sBarName -run $nMinRun -gap $nMaxGap -thr $nThreshold -less $nLess";
system($sCommand)== 0 or die "system '$sCommand' failed: $?";
my $sTilecoreOutput = get_tilecore_output();
if ($sTilecoreOutput){
   $sOutFile = $sTilecoreOutput unless($sOutFile);
   move($sTilecoreOutput, "$sCwd/$sOutFile") or die "Error: Could not move temporary tilecore output file: $!\n";
   chdir $sCwd;
}
else{
   my $sBedOut = join("", basename($sBarName, '.bar'), ".bed");
   chdir $sCwd;
   open BEDOUT, ">$sBedOut" or die "Error: couldn't write '$sBedOut': $!\n";
   close BEDOUT;
   warn "Warning: no output was found. This likely means no significant segments were detected, or there was a problem running tilecore\n";
}

#################
## SUBROUTINES ##
#################

# get_tilecore_output
#
# Gets the output file generated by tilecore
sub get_tilecore_output{
   opendir(DIR, ".");
   my @asFiles = grep(/\.bed$/,readdir(DIR));
   closedir(DIR);
   my $sReturn = @asFiles ? shift @asFiles : '';
   return $sReturn;
}


# INTERRUPT
#
# Interrupt routine, make sure we exit gracefully for tmp file cleanup
sub INTERRUPT{
   exit(1); # This will call END
}
