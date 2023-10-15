#!/usr/bin/env perl

# 21.09.2009 10:49:58 EDT
# Harm van Bakel <hvbakel@gmail.com>

# GLOBALS
$ENV{TMPDIR}    ||= "/tmp";     # location for tmp file storage
$ENV{SORT}      ||= "sort";     # the sort binary
$SIG{'INT'}=$SIG{'HUP'}=$SIG{'ABRT'}=$SIG{'QUIT'}=$SIG{'TRAP'}=$SIG{'STOP'}=\&INTERRUPT;

# MODULES
use strict;
use warnings;
use Getopt::Long;
use FileHandle;
use File::Temp qw(tempfile tempdir);
use ValidateFiletypes qw(check_gff);

# GET PARAMETERS
my $sHelp          = 0;
my $sArrayDataFile = '';
my $sProbeMapFile  = '';
GetOptions("help!"        => \$sHelp,
           "array-data:s" => \$sArrayDataFile,
           "probemap:s"   => \$sProbeMapFile);

# PRINT HELP
$sHelp = 1 unless($sArrayDataFile and $sProbeMapFile);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Script to map array data onto a genome

   Usage: $sScriptName -p <probemap> -a <array-data>
    
    -p <string>
      Probe mapping file in gff format. The identifier in this file
      must match the probe identifier used in the array data matrix
    -a <string>
      File containing the array data matrix. First column should contain
      the probe identifier and the first row should contain a unique
      description for each data colum that will be used as output file 
      names
    -help
      This help message
   
HELP
}


###########
## START ##
###########

# Check for proper gff formatting
my @asGffErrors = check_gff($sProbeMapFile);
if (@asGffErrors){
   unshift @asGffErrors, "Errors in gff file:";
   die join("\n", @asGffErrors), "\n";
}

# Read the probe mappings into a hash
my $rOligoPositions = get_oligo_positions($sProbeMapFile);

# map and write the sgr files
my ($nMapped, $nUnmapped) = map_sgr_files($sArrayDataFile, $rOligoPositions);
print STDERR "Mapped probes:   $nMapped\n";
print STDERR "Unmapped probes: $nUnmapped\n";



#################
## SUBROUTINES ##
#################

# get_oligo_positions($sGenomeDir)
#
# Get mappings for all oligo's
sub get_oligo_positions {
   my $sTrackFile = shift @_;
   my %hOligoMappings;
   
   # Read oligo mappings
   open TRACK, $sTrackFile or die "Error: could not open '$sTrackFile': $\n";
   while(<TRACK>){
      next if (/^\s*$/);
      next if (/^\s*#/);
      next if (/^\s*[Tt]rack/);
      s/\n//g; s/\r//g;
      s/\"//g; s/\'//g; # remove quotes!
      my ($sID, $sChr, $nStart, $nEnd) = (split /\t/)[8,0,3,4];
      my $nMid = $nStart + int( ($nStart-$nEnd+1) / 2 );
      push @{$hOligoMappings{uc($sID)}}, [$sChr, $nMid];
   }
   return(\%hOligoMappings);
}


# map_sgr_files($sInputFile, $rOligoPositions)
#
# Creates the outut graphs
sub map_sgr_files{
   my ($sInputFile, $rOligoPositions) = @_;
   my $nMapcount   = 0;
   my $nUnmapcount = 0;
   my %hResults;
   
   # Read the input file, create a hash entry for each result file
   my $fhInput = new FileHandle;
   $fhInput->open($sInputFile) or die "Error: could not '$sInputFile': $!\n";
   my $sHeader =  $fhInput->getline();
   $sHeader    =~ s/[\"\'\/\\\r\n]//g;
   $sHeader    =~ s/ /\-/g;
   my ($tmp, @asHeader) = split /\t/, $sHeader;
   my $nDataCols = scalar(@asHeader);
   
   # Check the header and create a temporary output file for each probe file
   my %hMapColToTmpfile;
   my %hMapNameToTempfile;
   for(my $i=0 ; $i< @asHeader ; $i++){
      my ($fhTmp, $sTmp) = tempfile('sgr-output-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1);
      die "Error: array data column names must be unique, but found duplicates for '$asHeader[$i]'\n" if (exists $hMapNameToTempfile{$asHeader[$i]});
      $hMapColToTmpfile{$i}              = $fhTmp;
      $hMapNameToTempfile{$asHeader[$i]} = $sTmp;
   }
   
   # Now read the file data and print the output to the temporary files
   while (my $sLine = $fhInput->getline()) {
      next if ($sLine =~ /^\s*$/);
      $sLine =~ s/[\n\r\"\']//g;
      my ($sID, @anData) = split /\t/, $sLine;
      die "Error: was expecting $nDataCols data columns on line $.\n" unless (scalar(@anData) == $nDataCols);
      if (exists $rOligoPositions->{$sID}){
         foreach my $rPosition (@{$rOligoPositions->{$sID}}){
            for (my $i=0 ; $i<@anData ; $i++){
               $hMapColToTmpfile{$i}->print(join("\t", @$rPosition, $anData[$i]), "\n");
            }
         }
         $nMapcount++;
      }
      else{
         $nUnmapcount++;
      }
   }
   $fhInput->close();
   
   # Finally, sort the output files and move them to the current dir
   foreach my $sOutputName (keys %hMapNameToTempfile){
      open OUT,  ">$sOutputName.sgr" or die "Error: can't write to file '$sOutputName': $!\n";
      open SORT, "$ENV{SORT} -S 2G -s -t '\t' -k1,1 -k2,2n $hMapNameToTempfile{$sOutputName} |" or die "Error: can't sort temporary file '$hMapNameToTempfile{$sOutputName}': $!\n";
      while (<SORT>){
         print OUT $_;
      }
      close SORT;
      close OUT;
   }
   
   return($nMapcount, $nUnmapcount);
}


# INTERRUPT
#
# Interrupt routine, make sure we exit gracefully for tmp file cleanup
sub INTERRUPT{
   exit(1); # This will call END
}
