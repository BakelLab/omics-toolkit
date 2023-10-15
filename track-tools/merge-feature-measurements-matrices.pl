#!/usr/bin/env perl

# Merge multiple feature measurement matrices into one big matrix

# MODULES
use strict;
use Getopt::Long;
use FileHandle;
use Cwd;


# GET ARGUMENTS
my $sHelp      = 0;
my $sIDonly    = 0;
my $sNames     = '';
my $sWeights   = '';
my $nSeparator = 0;
GetOptions("help!"       => \$sHelp,
           "id_only!"    => \$sIDonly,
           "names:s"     => \$sNames,
           "weights:s"   => \$sWeights,
           "sep!"        => \$nSeparator);


# PRINT HELP
$sHelp = 1 unless(@ARGV);
if ($sHelp) {
    die <<HELP

    merge-feature-measurements-matrices.pl [-names feature1,feature2 -weights weight1,weight2 ] matrix1 ... matrixN

    Arguments:
    -names <string>
      A comma-separated list of names that you want to attach to each matrix in the merged file.
      The list has to be in the same order as the specified matrix files
    -weights <string>
      A comma-separated list of weights between 0 and 1, one for each matrix
    -id_only
      Specify if you want to drop the strand and feature_type info in the merged file
      This makes things a bit easier for clustering
    -sep
      Include a separator (NA) between matrices
    -help
      This help message
      
HELP
}


###########
## START ##
###########

# Check arguments, parse names and check inputted matrix files
my @asNames;
my @anWeights;
&check_matrix_files(@ARGV);
if ($sNames){
   @asNames = &parse_names($sNames);
   die "Error: Number of names does not match the number of matrix files\n" unless (scalar(@asNames) == scalar(@ARGV));
}
if ($sWeights){
   @anWeights = &parse_weights($sWeights);
   die "Error: Number of weights does not match the number of matrix files\n" unless (scalar(@anWeights) == scalar(@ARGV));
}


# Start processing each matrix file
my %hAttributes;
my %hAttributesHeaders;
my %hIDs;
my %hData;
my %hDataHeaders;
foreach my $sFile (@ARGV){
   # Open the matrix file and process the header
   my $fhMatrix = new FileHandle;
   $fhMatrix->open($sFile) or die "Can't open the matrix file '$sFile': $!\n";
   
   # Process the file header; extract the feature attributes and the data header and keep both
   my $rAttributesHeader;
   my $sHeader           = $fhMatrix->getline();
   ($rAttributesHeader, $hDataHeaders{$sFile}) = &get_header_fields($sHeader, $sFile);
   my $nAttributesCols   = scalar(@$rAttributesHeader);
   my $nDataCols         = scalar(@{$hDataHeaders{$sFile}});
   foreach my $sName (@$rAttributesHeader){
      my $nID = scalar(keys(%hAttributesHeaders));
      $hAttributesHeaders{$sName} = $nID unless(exists($hAttributesHeaders{$sName}));
   }
   
   # Start processing the data
   while($_ = $fhMatrix->getline()){
      next if /^\s*#/;
      next if /^\s*$/;
      s/[\n\r]//g;
      my @anData = split /\t/, $_, -1;
      my $sID            = splice @anData, 0, 1;
      my (@asAttributes) = splice @anData, 0, $nAttributesCols;
      die "Error: Number of data points for feature '$sID' in '$sFile' does not match number of data headers\n" unless (scalar(@anData) == $nDataCols);
      
      # Keep track of feature IDs across files
      $hIDs{$sID}++;
      
      # Store all feature attributes; if a feature attribute was already set after processing a previous file,
      # make sure that the attributes are the same between files
      for (my $i=0 ; $i<$nAttributesCols ; $i++){
         if(exists($hAttributes{$sID}{$rAttributesHeader->[$i]})){
            unless ($hAttributes{$sID}{$rAttributesHeader->[$i]} eq $asAttributes[$i]){
               die "Error: attribute '$rAttributesHeader->[$i]' for feature '$sID' does not have the same value in all input files\n";
            }
         }
         else{
            $hAttributes{$sID}{$rAttributesHeader->[$i]} = $asAttributes[$i];
         }
      }
      
      # Store data values in following hash structure: matrixfile -> ID -> values
      # Make sure to check for uniqueness of IDs
      if (exists($hData{$sFile}{$sID})){
         die "Error: The feature '$sID' occurs more than once in matrix file '$sFile'\n";
      }
      else{
         $hData{$sFile}{$sID} = [@anData];
      }
   }
}


# Print the header of the output file
my @asSortedAttributes = sort {$hAttributesHeaders{$a} <=> $hAttributesHeaders{$b}} keys(%hAttributesHeaders);
my @asHeaderNames;
my @asHeaderWeights;
push @asHeaderNames, "YORF";
push @asHeaderWeights, "EWEIGHT" if (@anWeights);
unless ($sIDonly){
   push @asHeaderNames, @asSortedAttributes;
   if (@anWeights){
      foreach my $sAttributeName (@asSortedAttributes) {push @asHeaderWeights, '';};
   }
}
for(my $i=0 ; $i<@ARGV ; $i++){
   # Substitute filename with name if names array is filled
   my $sMatrixName = $ARGV[$i];
   if (@asNames){
      $sMatrixName = $asNames[$i];
   }
   else{
      # strip any extension and path from the name
      $sMatrixName =~ s/.*\///;
      $sMatrixName =~ s/\..*//;
   }

   foreach my $sHeaderField (@{$hDataHeaders{$ARGV[$i]}}){
      push @asHeaderNames, "$sMatrixName.$sHeaderField";
      push @asHeaderWeights, $anWeights[$i] if (@anWeights);
   }
   if ($nSeparator){
      push @asHeaderNames, 'NA';
      push @asHeaderWeights, 0 if (@anWeights);
   }
}
print join("\t", @asHeaderNames), "\n";
print join("\t", @asHeaderWeights), "\n" if (@anWeights);


# Make sure any EWEIGHT key is at the beginning of the list
my @asSortedFeatures = sort(keys(%hIDs));
if (exists $hIDs{'EWEIGHT'}){
   my $nPos = 0;
   for (my $i=0 ; $i<@asSortedFeatures ; $i++){
      if ($asSortedFeatures[$i] eq 'EWEIGHT'){
         $nPos = $i;
         last;
      }
   }
   if (@anWeights){
      splice(@asSortedFeatures, $nPos, 1);  # Specified weights trump weights already in the matrix files
   }
   else{
      unshift @asSortedFeatures, splice(@asSortedFeatures, $nPos, 1);
   }
}

# Print the feature attributes and data matrix
foreach my $sFeature (@asSortedFeatures){
   my @anDataLine;
   my $flMissingFeatures = 0;
   
   # Collect data values
   foreach my $sFile (@ARGV){
      # Check if combination of matrix and feature exist. If so, push to line-array to be printed, else add a set of NA's
      if (exists($hData{$sFile}{$sFeature})){
         push @anDataLine, @{$hData{$sFile}{$sFeature}};
      }
      else{
         foreach my $sMissing (@{$hDataHeaders{$sFile}}){ push @anDataLine, 'NA';};
      }
      push @anDataLine, 'NA' if ($nSeparator);
   }

   # Print data values with feature attributes
   if ($sIDonly){
      print join("\t", $sFeature, @anDataLine), "\n";
   }
   else{
      my @asAttributes;
      foreach my $sAttributeName (@asSortedAttributes){
         if(exists($hAttributes{$sFeature}{$sAttributeName})){
            push @asAttributes, $hAttributes{$sFeature}{$sAttributeName};
         }
         else{
            push @asAttributes, 'NA';
         }
      }
      print join("\t", $sFeature, @asAttributes, @anDataLine), "\n";
   }
}


#################
## SUBROUTINES ##
#################

# parse_names($sNames)
#
# Parse the names string and check if each name is unique
sub parse_names{
   my $sNames = shift @_;
   my @asReturn;
   my %hNames;
   
   # Split and check uniqueness
   my @asNames = split /,/, $sNames;
   foreach my $sName (@asNames){
      if (exists($hNames{$sName})){
         die "Error: the name '$sName' occurs more than once. All names must be unique\n";
      }
      else{
         push @asReturn, $sName;
         $hNames{$sName}++;
      }
   }
   return @asReturn;
}


# parse_weights($sWeights)
#
# Parse the weights string and checks if weights are numeric between 0 and 1
sub parse_weights{
   my $sWeights = shift @_;
   
   # Split and check values
   my @anWeights = split /,/, $sWeights;
   foreach my $sWeight (@anWeights){
      if ($sWeight =~ /\d+\.*\d*/){
         die "Weight must be a number between 0 and 1\n" unless ( ($sWeight>=0) and ($sWeight<=1) );
      }
      else{
         die "Weight must be a number between 0 and 1\n";
      }
   }
   return @anWeights;
}


# check_matrix_files(@asMatrixFiles)
#
# Remove any duplicate names from the matrix file input
sub check_matrix_files{
   my @asMatrixFiles = @_;
   
   # Go through each element and check uniqueness
   my %hFiles;
   foreach my $sFile (@asMatrixFiles){
      if (exists($hFiles{$sFile})){
         die "Error: the matrix file '$sFile' was specified more than once.\n";
      }
      else{
         $hFiles{$sFile}++;
      }
   }
   return 1;
}


# get_header_fields ($sHeader)
#
# Gets the constant and non-constant header part (i.e. the header of the data values itself)
sub get_header_fields {
   my ($sHeader, $sFile) = @_;
   my @asAttributesHeader;
   my @asDataHeader;
   
   # Parse the header
   $sHeader       =~ s/[\n\r]//g;
   my (@asHeader) =  split /\t/, $sHeader;
   if(lc($sHeader) =~ /^#*(geneid|yorf)\tfeature_type\tchr\tstrand\tparent_start\tparent_end\tparent_spliced_size\t/){
      @asAttributesHeader = @asHeader[1..6];
      @asDataHeader       = @asHeader[7..$#asHeader];
   }
   elsif(lc($sHeader) =~ /^#*(geneid|yorf)\tstrand\tfeature_type/){
      @asAttributesHeader = @asHeader[1..2];
      @asDataHeader       = @asHeader[3..$#asHeader];
   }
   elsif(lc($sHeader) =~ /^#*(geneid|yorf)\t/){
      @asDataHeader       = @asHeader[1..$#asHeader];
   }
   else{
      die "Error: matrix file '$sFile' does not have a valid header\n";
   }
   return (\@asAttributesHeader, \@asDataHeader);
}

