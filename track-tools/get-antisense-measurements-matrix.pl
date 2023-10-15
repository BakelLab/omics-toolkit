#!/usr/bin/env perl

# Get sense and antisense matrix files from forward and reverse strand matrices

# MODULES
use strict;
use Getopt::Long;
use FileHandle;
use Cwd;


# GET ARGUMENTS
my $sHelp     = 0;
my $sForward  = '';
my $sReverse  = '';
GetOptions("help!"       => \$sHelp,
           "forward:s"   => \$sForward,
           "reverse:s"   => \$sReverse);


# PRINT HELP
$sHelp = 1 unless($sForward and $sReverse);
if ($sHelp) {
    die <<HELP

    get-sense-measurements-matrix.pl 

    Arguments:
    -forward <string>
      The matrix file with the forward strand measurements
    -reverse <string>
      The matrix file with the reverse strand measurements
    -help
      This help message

HELP
}


###########
## START ##
###########

# Load sense and antisense features in a hash using geneID as key (must be unique otherwise throw error)
# Keep a separate hash with all gene identifiers in both files, we'll use this later to make sure that
# both supplied files contain all gene identifiers

my %hSense;
my %hAntisense;
my %hGeneIDs;

# Read the forward and reverse matrix files
my $sHeaderFwd = &read_matrix_file($sForward, 'fwd', \%hSense, \%hAntisense, \%hGeneIDs);
my $sHeaderRev = &read_matrix_file($sReverse, 'rev', \%hSense, \%hAntisense, \%hGeneIDs);

# Make sure that the forward and reverse strand matrix file headers match
die "Error: header mismatch for forward and reverse strand matrix files\n" unless ($sHeaderFwd eq $sHeaderRev);

# Now output the antisense matrix
print "$sHeaderFwd\n";
foreach my $sFeature (sort(keys(%hGeneIDs))){
   if ( exists $hAntisense{$sFeature} ){
      print join("\t", $sFeature, @{$hAntisense{$sFeature}}), "\n";
   }
   else{
      print "Warning: feature '$sFeature' is missing in one of the input files, feature skipped\n";
   }
}


#################
## SUBROUTINES ##
#################

# read_matrix_file ($sFile, $sType, \%hSense, \%hAntisense, \%hGeneIDs)
#
# Read the matrix file and put in sense and antisense hashes based on strand type
# Returns the file header
sub read_matrix_file{
   my ($sFile, $sType, $rhSense, $rhAntisense, $rhGeneIDs) = @_;

   # Open the matrix file and check if it has a valid header
   my $fhMatrix = new FileHandle;
   $fhMatrix->open($sFile) or die "Can't open the matrix file '$sFile': $!\n";
   my $sHeader  =  $fhMatrix->getline();
   $sHeader     =~ s/[\n\r]//g;
   $sHeader     =  lc($sHeader);
   die "Error: matrix file '$sFile' does not have a valid header\n" unless($sHeader =~ /^(#geneid|yorf)\tstrand\tfeature_type/);
   
   # Now start processing the data, split into sense and antisense matrices
   while($_ = $fhMatrix->getline()){
      next if /^\s*#/;
      next if /^\s*$/;
      s/[\n\r]//g;
      my ($sID, $sStrand, $sFeatureType, @anMatrix) = split /\t/, $_, -1;
      
      $rhGeneIDs->{$sID}++;
      if ($sType eq 'fwd'){
         if ($sStrand eq '+')   { $rhSense->{$sID}     = [($sStrand, $sFeatureType, @anMatrix)];   }
         elsif ($sStrand eq '.'){ $rhSense->{$sID}     = [($sStrand, $sFeatureType, @anMatrix)];   }
         elsif ($sStrand eq '-'){ $rhAntisense->{$sID} = [($sStrand, $sFeatureType, @anMatrix)];   }
         else                   { die "Unknown strand type '$sStrand' for feature '$sID'\n"; }
      }
      elsif ($sType eq 'rev'){
         if ($sStrand eq '+')   { $rhAntisense->{$sID} = [($sStrand, $sFeatureType, @anMatrix)];   }
         elsif ($sStrand eq '.'){ $rhAntisense->{$sID} = [($sStrand, $sFeatureType, @anMatrix)];   }
         elsif ($sStrand eq '-'){ $rhSense->{$sID}     = [($sStrand, $sFeatureType, @anMatrix)];   }
         else                   { die "Unknown strand type '$sStrand' for feature '$sID'\n"; }
      }
      else{
         die "Error: Unknown strand type '$sType' specified in subroutine 'read_matrix_file'\n";
      }
      
   }
   $fhMatrix->close();
   
   return $sHeader;
}
