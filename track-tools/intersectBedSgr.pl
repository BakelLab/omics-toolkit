#!/usr/bin/env perl

# 10.07.2011 16:19:02 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use ValidateFiletypes qw(check_bed);
use File::Temp qw(tempfile);

# GLOBALS
$ENV{BAR2TXT}      ||= "bar2txt";       # The bar2txt binary
$ENV{INTERSECTBED} ||= "intersectBed";  # The bedTools intersectBed binary
$ENV{TEMPDIR}      ||= '/tmp';          # Location for tmp file storage

# GET PARAMETERS
my $sHelp        = 0;
my $sBedFile  = "";
my $sSgrFile  = "";
my $sBarFile  = "";
my $flRescale = 0;
GetOptions("help!"    => \$sHelp,
           "bed:s"    => \$sBedFile,
           "sgr:s"    => \$sSgrFile,
           "bar:s"    => \$sBarFile,
           "rescale!" => \$flRescale);


# PRINT HELP
$sHelp = 1 unless($sBedFile and ($sSgrFile or $sBarFile));
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
   
   Get sgr or bar file data points that overlap bed features.
   
   Arguments:
    -bed <string>
      Bed file
    -sgr <string>
      Sgr file
    -bar <string>
      Bar file
    -rescale
      Rescale sgr coordinates relative to feature start
    -help
      This help message
   
HELP
}

##########
## MAIN ##
##########

# Check input
if ($sSgrFile and $sBarFile){
   die "Error: specify either an sgr or bar file as input, not both\n";
}

# Check bed file
my @asBedErrors = check_bed($sBedFile, 6);
if (@asBedErrors){
   unshift @asBedErrors, "The following errors were found in the input bed file:";
   die join("\n", @asBedErrors), "\n";
}

# Get a temporary sgr file to feed into intersectBed
my $sTmpSgr = $sSgrFile ? preprocess_sgr($sSgrFile) : preprocess_bar($sBarFile);

# Now run intersectBed

if ($flRescale){
   open INTERSECT, "$ENV{INTERSECTBED} -wo -a $sBedFile -b $sTmpSgr |" or die "Error: can't run intersectBed: $!\n";
   while (<INTERSECT>){
      next if (/^\s*$/);
      next if (/^ *#/);
      my @asLine = split(/\t/, $_, -1);
      my ($sChr,$nFeatureStart,$nFeatureEnd,$sFeatureID,$sFeatureStrand,$nSgrPos,$nSgrVal) = @asLine[0,1,2,3,5,7,9];
      my $nSgrPosNew = $sFeatureStrand eq '-' ? $nFeatureEnd - $nSgrPos : $nSgrPos - $nFeatureStart;
      print join ("\t", $sChr, $nSgrPosNew, $nSgrVal, $sFeatureID), "\n";
   }
   close INTERSECT;
}
else{
   open INTERSECT, "$ENV{INTERSECTBED} -u -a $sTmpSgr -b $sBedFile |" or die "Error: can't run intersectBed: $!\n";
   while (<INTERSECT>){
      next if (/^\s*$/);
      next if (/^ *#/);
      my @asLine = split(/\t/, $_, -1);
      print join("\t", @asLine[0,1,3]);
   }
   close INTERSECT;
}


#################
## SUBROUTINES ##
#################

# preprocess_bar
#
# preprocess bar file
sub preprocess_bar {
   my ($sBarFile) = @_;
   my ($fhTmpOut, $sTmpOut) = tempfile('intersectBedSgr-XXXXX', DIR=>$ENV{TEMPDIR}, UNLINK=>1);
   open BAR, "$ENV{BAR2TXT} -bar $sBarFile |" or die "Error: can't run bar2txt: $!\n";
   while (<BAR>){
      next if (/^ *#/);
      my (@asLine) = split("\t", $_, -1);
      print $fhTmpOut join("\t", $asLine[0], $asLine[1], $asLine[1]+1, $asLine[2]);
   }
   close BAR;
   $fhTmpOut->close();
   return $sTmpOut;
}

# preprocess_sgr
#
# preprocess sgr file
sub preprocess_sgr {
   my ($sSgrFile) = @_;
   my ($fhTmpOut, $sTmpOut) = tempfile('intersectBedSgr-XXXXX', DIR=>$ENV{TEMPDIR}, UNLINK=>1);
   open SGR, "$sSgrFile" or die "Error: can't open '$sSgrFile': $!\n";
   while (<SGR>){
      next if (/^\s*$/);
      next if (/^ *#/);
      my (@asLine) = split("\t", $_, -1);
      die "Error: inconsistent sgr format in '$sSgrFile', line $.\n" unless (@asLine==3 and $asLine[1]=~/^\d+$/);
      print $fhTmpOut join("\t", $asLine[0], $asLine[1], $asLine[1]+1, $asLine[2]);
   }
   close SGR;
   $fhTmpOut->close();
   return $sTmpOut;
}
