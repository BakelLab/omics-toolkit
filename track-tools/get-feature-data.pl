#!/usr/bin/env perl

# Intersect features with measurements based on genomic position mappings and
# get a binned matrix of measurements around (flanking) and inside features 
# that can be used for producing average plots or for clustering

# MODULES
use strict;
use Getopt::Long;
use FileHandle;
use Cwd;
use POSIX;
use File::Temp qw(tempfile tempdir);

# GLOBALS
$ENV{TMPDIR}         = "/sc/orga/scratch/vanbah01/tmp/";        # location for tmp file storage
$ENV{BEDTLS}         ||= "bedtools";    # location of bedtools binary
$ENV{BAR2TXT}        ||= "bar2txt";     # location of bar2txt binary
$ENV{SORT}           ||= "sort";        # Unix sort binary
$ENV{BIN_SIZE_RATIO} ||= 1;             # specifies lower limit of the feature or flanking
                                        # element size in relation to the bin size
                                        # i.e. specifying 10 means that with 3 bins, the 
                                        # feature or flanking size must be at least 30
$SIG{'INT'}=$SIG{'HUP'}=$SIG{'ABRT'}=$SIG{'QUIT'}=$SIG{'TRAP'}=$SIG{'STOP'}=\&INTERRUPT;

# GET ARGUMENTS
my $sHelp                = 0;
my $sFeaturesFile        = '';
my $sFeaturesFormat      = 'gff';
my $sDataFile            = '';
my $sDataFwd             = '';
my $sDataRev             = '';
my $sDataFormat          = 'bar';
my $sMappingType         = 'feature';
my $nBinCount            = 1;
my $sBinMeasure          = 'mean';
my $nFlankSize           = 500;
my $sExcludeFlankOverlap = 'none';
my $nFlankBufferSize     = 50;
my $nMinFeatureSize      = 0;
my $sOutputFormat        = 'relative';
my $sDataStrand          = 'both';
GetOptions("help!"                   => \$sHelp,
           "feature-file:s"          => \$sFeaturesFile,
           "feature-format:s"        => \$sFeaturesFormat,
           "data-file:s"             => \$sDataFile,
           "data-format:s"           => \$sDataFormat,
           "data-fwd:s"              => \$sDataFwd,
           "data-rev:s"              => \$sDataRev,
           "mapping-type:s"          => \$sMappingType,
           "bin-count:n"             => \$nBinCount,
           "bin-measure:s"           => \$sBinMeasure,
           "flank-size:n"            => \$nFlankSize,
	        "flank-buffer_size:n"     => \$nFlankBufferSize,
           "exclude-flank-overlap:s" => \$sExcludeFlankOverlap,
           "min-feature-size:n"      => \$nMinFeatureSize,
           "output-format:s"         => \$sOutputFormat,
           "strand:s"                => \$sDataStrand);


# PRINT HELP
$sHelp = 1 unless($sFeaturesFile and ($sDataFile or ($sDataFwd and $sDataRev)) and $nBinCount);
if ($sHelp) {
    my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
    die <<HELP

    $sScriptName [options] -feature-file <file> {-data-file <file> | -data-fwd <file> -data-rev <file>}

    Intersect genomic features with genome-wide measurements to plot data values relative to
    the feature coordinates. The resulting output can be used for cluster diagrams or making
    'averagogram' plots.

    Input:
    -feature-file <string>
      A tab-delimited file containing genome position information for the features
    -feature-format <string>
      Can be either gff or bed_basic. The bed_basic format only takes the first six
      fields into account (chr, start, end, ID, score, strand) and does not account 
      for exon structure!
      default: $sFeaturesFormat
    -data-file <string>
      The file containing the data values together with genomic positions
    -data-fwd <string>
      File containing data values for the forward/top/positive strand of the genome
    -data-rev <string>
      File containing data values for the reverse/bottom/negative strand of the genome
    -data-format <string>
      Can be either bar, sgr or bam. Bar files are binary mapping files produced by the
      affymetrix TAS software, sgr files are simple tab-delimited files with three columns:
      chromsome, position and data value. Bam files are NGS mapped read files in binary format.
      default: $sDataFormat
   
    Options:
    -strand <string>
      Specify strand to consider measurements from (only valid for bam files)
       both  :   Use measurements from both strands
       same  :   Require same strandedness as the feature
       anti  :   Require different strandedness
      default: $sDataStrand
    -mapping-type <string>
      Indicates how you want to use the annotation data for mapping, select one of:
       feature   :   mapping of complete features, excluding introns
       start     :   flanking only, relative to annotation start, e.g. TSS
       end       :   flanking only, relative to annotation end
       default: $sMappingType 
    -bin-count <integer>
      The number of regularly spaced bins to use for grouping the measurements in 
      the flanking regions and within features. This number must be at least 1
      default: $nBinCount
    -bin-measure <string>
      The type of measurement to take in each bin: cov, mean, median, sum or nfr (10th percentile)
      default: $sBinMeasure
    -flank-size <integer>
      The size of the upstream and downstream flanking regions to include in the output
      matrix. Set flank_size to zero to exclude flanking regions from the output matrix.
      default: $nFlankSize
    -flank-buffer-size <integer>
      Additional buffer to increase distance to neighboring features
      default: $nFlankBufferSize
    -exclude-flank-overlap
      Prevent flanking regions from overlapping neighboring features. This is for example
      useful for transcript mapping to exclude expression signals from neighboring features.
      Bins overlapping other features will be set to 'NA' to indicate missing data.
       both  :   Prevent overlap with features on either strand
       same  :   Prevent overlap with features on the same strand
       anti  :   Prevent overlap with features on the opposite strand
       none  :   Ignore overlapping features
      default: $sExcludeFlankOverlap
    -min-feature-size <integer>
      Only include features with a minimal size specified here in the output matrix. Set to
      zero to disable filtering on feature size
      default: $nMinFeatureSize
    -output-format <string>
      The output format, can be one of:
       id_only  :    Feature ID with data matrix
       basic    :    Feature ID, strand, type, with data matrix
       extended :    Feature ID, type, chr, strand, start, end, spliced-size, with data matrix
       relative :    Print relative coordinates to feature start or end (no bins)
    -help
      This help message
      
HELP
}


#################
## CHECK INPUT ##
#################

# Input files
die "Error: Can't find $sFeaturesFile at specified location\n"     unless (-e "$sFeaturesFile");

if ($sDataFile){
   die "Error: The -data-file option cannot be used in combination with -data-fwd or -data-rev\n" if ($sDataFwd or $sDataRev);
   die "Error: Can't find $sDataFile at specified location\n" unless (-e "$sDataFile");
}
elsif ($sDataFwd and $sDataRev){
   die "Error: Can't find $sDataFwd at specified location\n" unless (-e "$sDataFwd");
   die "Error: Can't find $sDataRev at specified location\n" unless (-e "$sDataRev");
}
else{
   die "Error: Missing require data files\n";
}

# Mapping type
my $flFlankOnly = 1;
$sMappingType = lc($sMappingType);
if ($sMappingType eq 'feature'){
   $flFlankOnly = 0;
}
elsif ($sMappingType eq 'start')  { 
   die "Error: Mapping type set to 'start' but no flank_size specified!\n" unless $nFlankSize;
}
elsif ($sMappingType eq 'end')    {
   die "Error: Mapping type set to 'end' but no flank_size specified!\n" unless $nFlankSize;
}
else { 
   die "Error: Unknown mapping type specified: $sMappingType\n"; 
}

# Other options
if ($nFlankSize){
   if ($nFlankSize< ($nBinCount*$ENV{BIN_SIZE_RATIO}) ){
      die join("\n", "Error: Too many bins specified relative to the size of the flanking region.",
                     " The size of the flanking region needs to be at least $ENV{BIN_SIZE_RATIO} times",
                     " greater than the number of bins to ensure sufficient measurements in each bin.");
   }
}
$sExcludeFlankOverlap = lc($sExcludeFlankOverlap);
die "Error: Unknown flank overlap type '$sExcludeFlankOverlap'\n" unless($sExcludeFlankOverlap =~ /^(both|same|anti|none)$/);
$sBinMeasure = lc($sBinMeasure);
die "Error: Unknown bin_measure type selected\n" unless ($sBinMeasure =~ /^(mean|median|sum|nfr|cov)$/);
$sOutputFormat = lc($sOutputFormat);
die "Error: Unknown output_format type selected\n" unless ($sOutputFormat =~ /^(id_only|basic|extended|relative)$/);
$sDataStrand = lc($sDataStrand);
die "Error: Unknown strand type selected\n" unless ($sDataStrand =~ /^(both|same|anti)$/);

# Make sure to only use coverage option when dealing with bam files
die "Error: Coverage option only available for bam files" if ( ($sBinMeasure eq 'cov') and ($sDataFormat ne 'bam') );

##########
## MAIN ##
##########

# Prepare the measurement data for the overlaps
print STDERR "Preparing temporary measurements file for overlaps\n";
my $sDataJoin = "";
$sDataFormat = lc($sDataFormat);
if ($sDataFormat eq 'bar'){
   $sDataJoin = $sDataFile ? &bar_to_join($sDataFile) : &bar_to_join($sDataFwd, $sDataRev);
}
elsif ($sDataFormat eq 'sgr'){
   $sDataJoin = $sDataFile ? &sgr_to_join($sDataFile) : &sgr_to_join($sDataFwd, $sDataRev);
}
elsif ($sDataFormat eq 'bam'){
   $sDataJoin = $sDataFile;
}
else{
   die "Error: Unknown measurement file format specified: $sDataFormat\n";
}


# Start processing the feature data
print STDERR "Preparing temporary features file for overlaps\n";
my $rFeaturesHash;
$sFeaturesFormat = lc($sFeaturesFormat);
if ($sFeaturesFormat eq 'gff'){
   $rFeaturesHash = &read_gff_features($sFeaturesFile, $flFlankOnly);
}
elsif ($sFeaturesFormat eq 'bed_basic'){
   $rFeaturesHash = &read_bed_basic_features($sFeaturesFile, $flFlankOnly);
}
elsif ($sFeaturesFormat eq 'bed'){
   die "Error: Bed format is not supported yet\n";
}
else{
   die "Error: Unknown feature file format specified: $sFeaturesFormat\n";
}

# Add start and end of closest flanking feature, if ExcludeOverlap is set in addition to a flanking region 
if ($nFlankSize and ($sExcludeFlankOverlap ne 'none')){
   &set_flank_positions($rFeaturesHash, $nFlankSize, $nFlankBufferSize, $sExcludeFlankOverlap);
}

# Now print the temporary features file for overlapping
my $sFeaturesJoin = &feature_hash_to_join($rFeaturesHash, $sMappingType, $nFlankSize);

# Run overlaps and retrieve the mappings in a hash table
print STDERR "Running join\n";
my $sJoinOutput = &get_join_mapping($sFeaturesJoin, $sDataJoin, $nFlankSize, $sDataFormat, $sDataStrand, $sBinMeasure);

# Prepare the output matrix
if ($sOutputFormat  ne "relative"){
   print STDERR "Printing output matrix\n";
   &print_matrix($sJoinOutput, $rFeaturesHash, $nBinCount, $nFlankSize, $flFlankOnly, $sBinMeasure, $sOutputFormat);
}
else{
   print STDERR "Printing relative coordinates\n";
   &print_relative($sJoinOutput, $rFeaturesHash, $sMappingType,$nFlankSize);
}



#################
## SUBROUTINES ##
#################

# print_features_hash($rFeaturesHash)
#
# Print the complete feature hash structure for debugging purposes
sub print_features_hash {
   my $rFeaturesHash = shift @_;
   foreach my $sFeature (keys(%$rFeaturesHash)){
      foreach my $sKey (sort(keys( %{$rFeaturesHash->{$sFeature}}))) {
         if ($sKey eq 'features'){
            foreach my $rPair (@{$rFeaturesHash->{$sFeature}{$sKey}}) {
               print join ("\t", $sFeature, $sKey, join('-', @$rPair)), "\n";
            }
         }
         else{
            print join("\t", $sFeature, $sKey, $rFeaturesHash->{$sFeature}{$sKey}), "\n";
         }
      }
   }
}


# bar_to_join($sBarFile)
#
# Convert a bar file to a file that is suitable for mapping with join
# Returns the name of the temporary join file
sub bar_to_join {
   my ($sBarFwd, $sBarRev) = @_;
   
   # Make temporary output file
   my ($fhTmp, $sTmp) = tempfile('getmatrix-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1);
   
   # Process forward strand data file
   die "Error: expected binary file for 'bar' measurement format. Please check -data-file or -data-fwd\n" unless (-B $sBarFwd);
   open CONVERT, "$ENV{BAR2TXT} -bar $sBarFwd |" or die "Error: Can't convert the '$sBarFwd': $!\n";
   while (<CONVERT>){
      s/[\n\r]//g;
      if (/^chr(\d+|MT|mt|Mt|Mito|mito|X|x)/){                       # Match only 'real' chromosomes, skip controls
         my ($sChr, $sPos, $sVal) = split /\t/;
         print $fhTmp join("\t", $sChr, $sPos, $sPos+1, '', $sVal,'+'), "\n";
      }
   }
   close CONVERT;
   
   # Process reverse strand data file
   if ($sBarRev){
      die "Error: expected binary file for 'bar' measurement format. Please check -data-rev\n" unless (-B $sBarRev);
      open CONVERT, "$ENV{BAR2TXT} -bar $sBarRev |" or die "Error: Can't convert '$sBarRev': $!\n";
      while (<CONVERT>){
         s/[\n\r]//g;
         if (/^chr(\d+|MT|mt|Mt|Mito|mito|X|x)/){                       # Match only 'real' chromosomes, skip controls
            my ($sChr, $sPos, $sVal) = split /\t/;
            print $fhTmp join("\t", $sChr, $sPos, $sPos+1, '', $sVal, '-'), "\n";
         }
      }
      close CONVERT;
   }
   close $fhTmp;
   return($sTmp);
}


# sgr_to_join($sGprFile)
#
# Convert a sgr file to a file that is suitable for mapping with join
# Returns the name of the temporary join file
sub sgr_to_join{
   my ($sSgrFwd, $sSgrRev) = @_;
   
   # Make temporary output file
   my ($fhTmp, $sTmp) = tempfile('getmatrix-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1);
   
   # Process forward strand data file
   my $nCount = 1;
   die "Error: expected text file for 'sgr' measurement format. Please check -data-file or -data-fwd\n" if (-B $sSgrFwd);
   open CONVERT, "$sSgrFwd" or die "Error: Can't open '$sSgrFwd' for processing: $!\n";
   while (<CONVERT>){
      next if /^\s*$/;
      next if /^\s*#/;
      s/[\n\r]//g;
      my (@asLine) = split /\t/;
      if (@asLine >= 3){
         my ($sChr, $sPos, $sVal) = @asLine[0,1,2];
         die "Error: position field must be an integer in '$sSgrFwd' line $.\n" unless ($sPos =~ /^\d+$/);
         print $fhTmp join("\t", $sChr, $sPos, $sPos+1, "F" .$nCount++, $sVal, '+'), "\n";
      }
      else{
         die "Error: Insufficient fields for an sgr file: '$sSgrRev' line $.\n";
      }
   }
   close CONVERT;
   
   # Process reverse strand data file
   $nCount = 1;
   if ($sSgrRev){
      die "Error: expected text file for 'sgr' measurement format. Please check -data-rev\n" if (-B $sSgrRev);
      open CONVERT, "$sSgrRev" or die "Error: Can't open '$sSgrRev' for processing: $!\n";
      while (<CONVERT>){
         next if /^\s*$/;
         next if /^\s*#/;
         s/[\n\r]//g;
         my (@asLine) = split /\t/;
         if (@asLine >= 3){
            my ($sChr, $sPos, $sVal) = @asLine[0,1,2];
            die "Error: position field must be an integer in '$sSgrRev' line $.\n" unless ($sPos =~ /^\d+$/);
            print $fhTmp join("\t", $sChr, $sPos, $sPos+1, "R". $nCount++, $sVal, '-'), "\n";
         }
         else{
            die "Error: Insufficient fields for an sgr file: '$sSgrRev' line $.\n";
         }
      }
      close CONVERT;
   }
   
   close $fhTmp;
   return($sTmp);  
}


# read_gff_features($sFeatureFile)
#
# Reads the gff features into a hash
# returns hash with:  -> chromosome          : Chromosome
#                     -> strand              : Strand (+ or - or .)
#                     -> parent_start        : Start of first subfeature
#                     -> parent_end          : End of last subfeature
#                     -> parent_spliced_size : size of complete annotation without introns
#                     -> features            : array of arrays with start, end (only if FlankOnly = 0);
sub read_gff_features{
   my ($sFeaturesFile, $flFlankOnly) = @_;
   my %hFeatures;
   
   # Check GFF header, we'll support gff1 and gff2 for now
   my $nGffVersion = 1;
   my $fhGff       = new FileHandle;
   $fhGff->open($sFeaturesFile) or die "Can't open $sFeaturesFile: $!\n";
   my $sHeader     = $fhGff->getline();
   $sHeader        =~ s/#//g;
   if ($sHeader =~ /gff-version\t(1|2)/) { $nGffVersion = $1; }
   else  { 
      die join("\n", "Error: Missing gff header or unsupported gff version. Only gff versions 1 and 2 are currently supported",
                     " Example of gff header:", " ##gff-version <tab> 2"), "\n";
   }
   
   # Start reading the annotations in the gff file
   while($_ = $fhGff->getline()){
      next if /^\s*$/;
      next if /^track\s/;   # We need to skip the track lines in a ucsc/igb formatted gff file
      s/[\n\r]//g;
      my ($sChr, $sSrc, $sType, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sID) = split /\t/;
      $nStart--; # Make start zero-based (as in bed files, makes calculations easier!)
      
      # Make sure that we can extract an identifier from each line and do some other basic checks on gff parsing
      $sID = &extract_gff2_id($sID) if ($nGffVersion == 2);
      die "Error: line $. in '$sFeaturesFile' is missing an identifier (expecting GFF format)\n" unless $sID;
      die "Error: unknown strand '$sStrand' for feature '$sID', line $. in '$sFeaturesFile'\n" unless ($sStrand =~ /[+-.]/);
      die "Error: expected number for start position on line $. in '$sFeaturesFile'\n" unless ($nStart =~ /^\d+$/);
      die "Error: expected number for end position on line $. in '$sFeaturesFile'\n" unless ($nEnd =~ /^\d+$/);
      
      # Load the hash, note that we'll skip the array of features if flanking_only is set
      ($nStart, $nEnd) = ($nEnd, $nStart) if ($nStart > $nEnd); # Make sure that the start is always smaller than the end
      if (exists($hFeatures{$sID})) {
         die "Error: feature '$sID' maps to multiple chromosomes in the gff file\n"        unless ($hFeatures{$sID}{chr} eq $sChr);
         die "Error: feature '$sID' has mappings to both the forward and reverse strand\n" unless ($hFeatures{$sID}{strand} eq $sStrand);   
         $hFeatures{$sID}{parent_start} = $nStart if ($nStart < $hFeatures{$sID}{parent_start});
         $hFeatures{$sID}{parent_end}   = $nEnd if ($nEnd > $hFeatures{$sID}{parent_end});
         $hFeatures{$sID}{parent_spliced_size} += $nEnd - $nStart + 1;
         push @{$hFeatures{$sID}{features}}, [($nStart, $nEnd)] unless ($flFlankOnly);
      }
      else{
         $hFeatures{$sID}{chr}                 = $sChr;
         $hFeatures{$sID}{strand}              = $sStrand;
         $hFeatures{$sID}{type}                = $sType;
         $hFeatures{$sID}{parent_start}        = $nStart;
         $hFeatures{$sID}{parent_end}          = $nEnd;
         $hFeatures{$sID}{parent_spliced_size} = $nEnd - $nStart;
         push @{$hFeatures{$sID}{features}}, [($nStart, $nEnd)] unless ($flFlankOnly);
      }
   }
   $fhGff->close();
   return(\%hFeatures);
}


# extract_gff2_id($sID)
#
# Extract the gene_id from a gff2 identifier
# Returns the identifier or an empty string if the gene_id could not be extracted
sub extract_gff2_id{
   my $sFields  = shift @_;
   my @asFields = split /\s*;\s*/, $sFields;
   my $sGeneID = '';
   foreach my $sField (@asFields){
      $sField =~ s/["']//g;
      my ($sKey, @asRest) = split /\s+/, $sField;
      $sGeneID = join('_', @asRest) if (lc($sKey) =~ 'gene');
   }
   return $sGeneID;
}


# read_bed_basic_features($sFeatureFile)
#
# Reads whole bed features (only looking at feature start and end) into a hash
# returns hash with:  -> chromosome          : Chromosome
#                     -> strand              : Strand (+ or - or .)
#                     -> parent_start        : Start of first subfeature
#                     -> parent_end          : End of last subfeature
#                     -> parent_spliced_size : size of complete annotation without introns
#                     -> features            : array of arrays with start, end (only if FlankOnly = 0);
sub read_bed_basic_features{
   my ($sFeatureFile, $flFlankOnly) = @_;
   my %hFeatures;
   
   # Start reading the annotations in the Bed file
   my $nCount = 1;
   my $fhBed       = new FileHandle;
   $fhBed->open($sFeatureFile) or die "Can't open $sFeatureFile: $!\n";  
   while($_ = $fhBed->getline()){
      next if /^\s*$/;
      next if /^track\s/;   # We need to skip the track lines in a ucsc/igb formatted bed file
      s/[\n\r]//g;
      my (@asLine) = split /\t/;
      die "Error: not enough fields for bed_basic format on line $.\n" unless(@asLine >=4);
      my ($sChr, $nStart, $nEnd, $sID, $nScore, $sStrand, @asRest) = @asLine;
      
      # Make sure that we can extract an identifier from each line and do some other basic checks on bed parsing
      die "Error: line $nCount in '$sFeaturesFile' is missing an identifier (expecting bed_basic format)\n" unless $sID;
      die "Error: expected number for start position on line $nCount in '$sFeaturesFile'\n" unless ($nStart =~ /\d+/);
      die "Error: expected number for end position on line $nCount in '$sFeaturesFile'\n" unless ($nEnd =~ /\d+/);
      
      # Load the hash, note that we'll skip the array of features if flanking_only is set
      ($nStart, $nEnd) = ($nEnd, $nStart) if ($nStart > $nEnd); # Make sure that the start is always smaller than the end
      if (exists($hFeatures{$sID})) {
         die "Error: feature '$sID' maps to multiple chromosomes in the gff file\n"        unless ($hFeatures{$sID}{chr} eq $sChr);
         $hFeatures{$sID}{parent_start} = $nStart if ($nStart < $hFeatures{$sID}{parent_start});
         $hFeatures{$sID}{parent_end}   = $nEnd   if ($nEnd   > $hFeatures{$sID}{parent_end});
         $hFeatures{$sID}{parent_spliced_size} += $nEnd - $nStart;
         push @{$hFeatures{$sID}{features}}, [($nStart, $nEnd)] unless ($flFlankOnly);
      }
      else{
         $hFeatures{$sID}{chr}                 = $sChr;
         $hFeatures{$sID}{strand}              = $sStrand ? $sStrand : '+';
         $hFeatures{$sID}{type}                = 'gene';
         $hFeatures{$sID}{parent_start}        = $nStart;
         $hFeatures{$sID}{parent_end}          = $nEnd;
         $hFeatures{$sID}{parent_spliced_size} = $nEnd - $nStart;
         push @{$hFeatures{$sID}{features}}, [($nStart, $nEnd)] unless ($flFlankOnly);
      }
      $nCount++;
   }
   $fhBed->close();
   return(\%hFeatures);
}



# set_flank_positions($rhFeatures, $nFlankSize, $nFlankBufferSize, $nConsiderStrands)
#
# Set the start and end positions for the flanking segments, excluding overlap with neigboring features
sub set_flank_positions{
   my ($rhFeatures, $nFlankSize, $nFlankBufferSize, $sExcludeFlankStrand) = @_;
   
   # Prepare a new array for processing, needs to have id, chromosome, strand, feature_start, feature_end
   # Sort the array descending on everything
   my @aaPositions;
   foreach my $sFt (keys(%$rhFeatures)){
      push @aaPositions, [($sFt, $rhFeatures->{$sFt}{chr}, $rhFeatures->{$sFt}{strand}, 
                           $rhFeatures->{$sFt}{parent_start}, $rhFeatures->{$sFt}{parent_end})]
   }
   
   # Find upstream flanking position
   # Sort the array by chromosome and start position for each feature, loop through all the features before 
   # it in the list on the same chromosome and pick the feature with the highest end as flank_start
   @aaPositions = sort { lc $a->[1] cmp lc $b->[1] || $a->[3] <=> $b->[3] } @aaPositions;
   for (my $i=0 ; $i<@aaPositions ; $i++){
      my ($sFt, $sChr, $sStr, $nStart, $nEnd) = @{$aaPositions[$i]};
      my $nFlankStart = 1;
      if ($i>0){
         for (my $j=$i-1; $j>=0 ; $j--){
            last if (  $aaPositions[$j][1] ne $sChr );                                           # only consider the same chromosome            
            next if ( ($aaPositions[$j][2] ne $sStr) and ($sExcludeFlankStrand eq 'same') );     # only consider same strand, if specified
            next if ( ($aaPositions[$j][2] eq $sStr) and ($sExcludeFlankStrand eq 'anti') );     # only consider opposite strand, if specified
            $nFlankStart = $aaPositions[$j][4]+$nFlankBufferSize  if ($aaPositions[$j][4]+$nFlankBufferSize>$nFlankStart);
         }
      }
      $nFlankStart = $nStart if ($nFlankStart>$nStart);                             # for overlapping features
      $rhFeatures->{$sFt}{flank_start} = $nFlankStart;
   }
   
   # Find downstream flanking position
   # Sort the array by chromosome and end position for each feature, loop through all the features after it in 
   # the list on the same chromosome and pick the feature with the smallest start as flank_end.
   @aaPositions = sort { lc $a->[1] cmp lc $b->[1] || $a->[4] <=> $b->[4] } @aaPositions;
   for (my $i=0 ; $i<@aaPositions ; $i++){
      my ($sFt, $sChr, $sStr, $nStart, $nEnd) = @{$aaPositions[$i]};
      my $nFlankEnd = 999999999999;                                                              # Pick a number > the biggest chromosome size possible
      if ($i<$#aaPositions){
         for (my $j=$i+1 ; $j<@aaPositions ; $j++){
            last if (  $aaPositions[$j][1] ne $sChr );                                           # only consider the same chromosome
            next if ( ($aaPositions[$j][2] ne $sStr) and ($sExcludeFlankStrand eq 'same'));      # only consider same strand, if specified
            next if ( ($aaPositions[$j][2] eq $sStr) and ($sExcludeFlankStrand eq 'anti'));      # only consider opposite strand, if specified
            $nFlankEnd = $aaPositions[$j][3]-$nFlankBufferSize  if ($aaPositions[$j][3]-$nFlankBufferSize<$nFlankEnd);
         }
      }
      $nFlankEnd = $nEnd if ($nFlankEnd<$nEnd);                             # for overlapping features
      $rhFeatures->{$sFt}{flank_end} = $nFlankEnd;
   }  
}


# feature_hash_to_join($rhFeatures)
#
# Now we use the feature hash to produce all elements for the file that we will map through join
# Returns the name of the temporary join file
sub feature_hash_to_join{
   my ($hFeatures, $sMappingType, $nFlankSize) = @_;
   
   # Open output file and loop through all features
   my ($fhTmp, $sTmp) = tempfile('getmatrix-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1);
   foreach my $sFt (keys(%$hFeatures)){
      my ($sChr, $sStr, $sType, $nParentStart, $nParentEnd, $nParentSize, $nFlankStart, $nFlankEnd) = 
         @{$hFeatures->{$sFt}}{('chr','strand','type','parent_start','parent_end','parent_spliced_size','flank_start','flank_end')};   # assign locals for readability
      
      # Print the flanking regions (if requested)
      if ($nFlankSize){
         my ($nLowerStart, $nLowerEnd, $nUpperStart, $nUpperEnd, $sLowerType, $sUpperType) = (0,0,0,0,'','');
         
         # Get the start and end coordinates of the lower and upper flanking regions
         if ($sMappingType eq 'start'){
            if ($sStr eq '-'){
               $nLowerStart = $nParentEnd - $nFlankSize;
               $nLowerEnd   = $nParentEnd;
               $nUpperStart = $nParentEnd;
               $nUpperEnd   = $nParentEnd + $nFlankSize;
            }
            else{
               $nLowerStart = $nParentStart - $nFlankSize;
               $nLowerEnd   = $nParentStart;
               $nUpperStart = $nParentStart;
               $nUpperEnd   = $nParentStart+$nFlankSize;
            }
         }
         elsif ($sMappingType eq 'end'){
            if ($sStr eq '-'){
               $nLowerStart = $nParentStart - $nFlankSize;
               $nLowerEnd   = $nParentStart;
               $nUpperStart = $nParentStart;
               $nUpperEnd   = $nParentStart + $nFlankSize;
            }
            else{
               $nLowerStart = $nParentEnd - $nFlankSize;
               $nLowerEnd   = $nParentEnd;
               $nUpperStart = $nParentEnd;
               $nUpperEnd   = $nParentEnd + $nFlankSize;
            }
         }
         else{
               $nLowerStart = $nParentStart - $nFlankSize;
               $nLowerEnd   = $nParentStart;
               $nUpperStart = $nParentEnd;
               $nUpperEnd   = $nParentEnd + $nFlankSize;
         }
         
         # For the negative strand the upper flanking region is upstream the feature and the lower flank is downstream
         # For the positive strand the order is reversed
         if ($sStr eq '-'){
            $sLowerType = 'd';
            $sUpperType = 'u';
         }
         else{
            $sLowerType = 'u';
            $sUpperType = 'd';
         }
         
         # Correct starts and ends for flanking regions if they were defined
         $nLowerStart = 1 if ($nLowerStart<1);
         $nLowerStart = $nFlankStart if ($nFlankStart and ($nFlankStart>$nLowerStart) );
         $nUpperEnd   = $nFlankEnd if ($nFlankEnd and ($nFlankEnd<$nUpperEnd) );
         
         # Now acutally print the flanking regions
         print $fhTmp join("\t", $sChr, $nLowerStart, $nLowerEnd, $sFt, $sLowerType, $sStr, 0), "\n";
         print $fhTmp join("\t", $sChr, $nUpperStart, $nUpperEnd, $sFt, $sUpperType, $sStr, 0), "\n";
      }
      
      # Print the features
      if ($sMappingType eq 'feature'){
         # Write all the segments, for + strand, sort by start, ascending and print each one with right offset
         # For - strand features, sort by end, descending and print each one with the right offset
         my @aaSegments = @{$hFeatures->{$sFt}{features}};
         if ($sStr eq '-') { @aaSegments = sort { $b->[2] <=> $a->[2] | $b->[1] <=> $a->[1] } @aaSegments }
         else              { @aaSegments = sort { $a->[1] <=> $b->[1] | $a->[2] <=> $b->[2] } @aaSegments }
         
         # Now print all the segments with the correct offset to the previous element!
         my $nOffset = 0;
         foreach my $rSegment (@aaSegments){
            print $fhTmp join("\t", $sChr, $rSegment->[0], $rSegment->[1], $sFt, 'f', $sStr, $nOffset), "\n";
            $nOffset += $rSegment->[1] - $rSegment->[0];
         }
      }
   }
   close $fhTmp;
   return $sTmp;
}


# get_join_mapping
#
# Now use the two temporary join mapping files and run bedtools to get overlaps
# We'll parse the output as it comes out and put it in a new hash
# feature_name -> join_type -> [relative_position, value]
sub get_join_mapping{
   my ($sFeatures, $sData, $nFlankSize, $sFormat, $sOvlStrand, $sType) = @_;
   my ($fhTmp, $sTmp) = tempfile('getmatrix-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1);
   
   # Set strandedness for overlaps
   my $sStrandedness = '';
   $sStrandedness    = '-s' if ($sOvlStrand eq 'same');
   $sStrandedness    = '-S' if ($sOvlStrand eq 'anti');
   
   # Run overlaps
   if ($sFormat eq 'bam'){
      # Run coverageBed for bam input
      open BEDTOOLS, "$ENV{BEDTLS} coverage $sStrandedness -split -d -abam $sData -b $sFeatures |" or die "Error: can't run join: $!\n";
      while (<BEDTOOLS>){
         next if /^\s*$/;
         s/[\n\r]//g;
         my ($sFchr, $nFstart, $nFend, $sFt, $sJoinType, $sFstr, $nFoffset, $nPoffset, $nPint) = split /\t/;
         my $nRelPos;
         if ($sFstr eq '-') { 
            if ($sJoinType eq 'u'){ $nRelPos = $nFlankSize - ($nPoffset);                     }  # for upstream, the position must be relative to the specified flank start!!
            else                  { $nRelPos = ($nFend - $nFstart) - $nPoffset + $nFoffset;   }
         }
         else { 
            if ($sJoinType eq 'u'){ $nRelPos = ($nFend - $nFstart) - $nFlankSize + $nPoffset; }  # for upstream, the position must be relative to the specified flank start!!
            else                  { $nRelPos = $nPoffset + $nFoffset;                         }
         }
         print $fhTmp join(";", $sFt, $sJoinType, $nRelPos, $nPint), "\n" unless (($nPint==0) and ($sType eq 'cov'));
      }
      close BEDTOOLS;
   }
   else{
      # Run intersectBed for all other input formats
      open BEDTOOLS, "$ENV{BEDTLS} intersect $sStrandedness -wo -a $sFeatures -b $sData |" or die "Error: can't run join: $!\n";
      while (<BEDTOOLS>){
         next if /^\s*$/;
         s/[\n\r]//g;
         my ($sFchr, $nFstart, $nFend, $sFt, $sJoinType, $sFstr, $nOffset, $sPchr, $nPstart, $nPend, $nPname, $nPint, $sPstr) = split /\t/;
         my $nRelPos;
         if ($sFstr eq '-') { 
            if ($sJoinType eq 'u'){ $nRelPos = $nFlankSize - ($nPstart-$nFstart); }  # for upstream, the position must be relative to the specified flank start!!
            else                  { $nRelPos = $nFend   - $nPstart + $nOffset;    }
         }
         else { 
            if ($sJoinType eq 'u'){ $nRelPos = $nFlankSize - ($nFend - $nPstart); }  # for upstream, the position must be relative to the specified flank start!!
            else                  { $nRelPos = $nPstart - $nFstart + $nOffset;    }
         }
         print $fhTmp join(";", $sFt, $sJoinType, $nRelPos, $nPint), "\n";
      }
      close BEDTOOLS;
   }
   
   # Close file and return name
   close $fhTmp;
   return($sTmp);
}


# print_matrix
#
# Now use the join mappings to generate the data matrix
sub print_matrix{
   my ($sJoin, $rFeatures, $nBinCount, $nFlankSize, $flFlankOnly, $sBinMeasure, $sOutputFormat) = @_;
   my %hPrintedFeatures;
   
   # Print a nice header
   my @asSegmentTypes;
   push @asSegmentTypes, 'up' if ($nFlankSize);
   push @asSegmentTypes, 'ft'  unless ($flFlankOnly);
   push @asSegmentTypes, 'dn' if ($nFlankSize);
   my @asHeader  = ('#geneid');
   push @asHeader, ('strand', 'feature_type') if ($sOutputFormat eq 'basic');
   push @asHeader, ('feature_type', 'chr', 'strand', 'parent_start', 'parent_end', 'parent_spliced_size') if ($sOutputFormat eq 'extended');
   my $nIntervalSize = $nFlankSize / $nBinCount;
   foreach my $sSegmentType (@asSegmentTypes){
      for (my $i=0 ; $i<$nBinCount ; $i++){
         if ($sSegmentType eq 'up'){
            push @asHeader, ( ($i*$nIntervalSize) + int($nIntervalSize/2) ) - $nFlankSize;
         }
         elsif ($sSegmentType eq 'dn'){
            push @asHeader, ( ($i*$nIntervalSize) + int($nIntervalSize/2) );
         }
         else{
            push @asHeader, "$sSegmentType.$i";
         }
      }
   }
   print join("\t", @asHeader), "\n";
   
   # Print the data
   my $sLastFeature = '';
   my %hFeatureData = {u=>{}, d=>{}, f=>{}};
   open JOIN, "$ENV{SORT} -t ';' -k1,1 $sJoin |" or die "Error: can't open join output file: $!\n";
   while (<JOIN>){
      s/[\n\r]//g;
      my ($sFt, $sJoinType, $nPos, $nVal) = split /;/;
      $sLastFeature = $sFt unless ($sLastFeature);
      if ($sFt eq $sLastFeature){
         $hFeatureData{$sJoinType}{$nPos} = $nVal;
      }
      else{
         # Print matrix for previous feature
         print &get_matrix_line($sLastFeature, \%hFeatureData, $rFeatures, $nBinCount, $nFlankSize, $sBinMeasure, $sOutputFormat) . "\n";
         $hPrintedFeatures{$sLastFeature}++;   # Keep the IDs of all features that we're printing
         
         # Prepare for next feature
         $sLastFeature  = $sFt;
         %hFeatureData  = {u=>{}, d=>{}, f=>{}};     # Clear the hash table
         $hFeatureData{$sJoinType}{$nPos} = $nVal;
      }
   }
   # We still need to deal with the last feature here
   print &get_matrix_line($sLastFeature, \%hFeatureData, $rFeatures, $nBinCount, $nFlankSize, $sBinMeasure, $sOutputFormat) . "\n";
   $hPrintedFeatures{$sLastFeature}++;   # Keep the IDs of all features that we're printing
   close JOIN;
   
   # Finally, report which features were skipped during printing
   my @asSkippedFeatures;
   foreach my $sFt (sort(keys(%$rFeatures))){ push @asSkippedFeatures, $sFt unless(exists($hPrintedFeatures{$sFt})); }
   if (@asSkippedFeatures){
      print STDERR join(' ', 'Warning: The following features were skipped in the output matrix: ', join(",", @asSkippedFeatures)), "\n";
   }
}


# get_matrix_line
#
# Returns a single formatted line for the measurement matrix
sub get_matrix_line{
   my ($sFeature, $rFeatureData, $rFeatures, $nBinCount, $nFlankSize, $sBinMeasure, $sOutputFormat) = @_;
   
   my @anAveragedBins;
   push @anAveragedBins, &get_measurement_bins($rFeatureData->{u}, $nBinCount, $nFlankSize, $sBinMeasure) if ($nFlankSize);
   push @anAveragedBins, &get_measurement_bins($rFeatureData->{f}, $nBinCount, $rFeatures->{$sFeature}{parent_spliced_size}, $sBinMeasure) unless ($flFlankOnly);
   push @anAveragedBins, &get_measurement_bins($rFeatureData->{d}, $nBinCount, $nFlankSize, $sBinMeasure) if ($nFlankSize);
   
   if ($sOutputFormat eq 'id_only'){
      return join("\t", $sFeature, @anAveragedBins);
   }
   elsif ($sOutputFormat eq 'basic'){
      return join("\t", $sFeature, $rFeatures->{$sFeature}{strand}, $rFeatures->{$sFeature}{type}, @anAveragedBins);
   }
   else{
      return join("\t", $sFeature, $rFeatures->{$sFeature}{'type'}, $rFeatures->{$sFeature}{'chr'}, $rFeatures->{$sFeature}{'strand'}, $rFeatures->{$sFeature}{'parent_start'}, $rFeatures->{$sFeature}{'parent_end'}, $rFeatures->{$sFeature}{'parent_spliced_size'}, @anAveragedBins);
   }
}


# get_measurement_bins
#
# Takes a hash of probe positions and values and bins it according to the bin size
sub get_measurement_bins{
   my ($rPositions, $nBinCount, $nSegmentSize, $sType) = @_;
   my %hBins;
   my @asResult;

   my $nBinSize = $nSegmentSize/($nBinCount);

   # Assign all probe values to intervals
   foreach my $nPosition (sort {$a <=> $b} keys(%$rPositions) ){
MAP:  for (my $i=0 ; $i<=$nBinCount ; $i++ ){
         if (($i*$nBinSize)>$nPosition ){
            push @{$hBins{$i-1}}, $rPositions->{$nPosition};
            last MAP;
         }
      }
   }

   # Now go through all the intervals and print mean values or NA if not present
   for (my $i=0 ; $i<$nBinCount ; $i++){
      if (exists($hBins{$i})){
         my $nBinMeasure = 0;
         if    ($sType eq 'mean')   { $nBinMeasure = &mean(@{$hBins{$i}})          }
         elsif ($sType eq 'median') { $nBinMeasure = &median(@{$hBins{$i}})        }
         elsif ($sType eq 'sum')    { $nBinMeasure = &sum(@{$hBins{$i}})           }
         elsif ($sType eq 'cov')    { $nBinMeasure = &sum(@{$hBins{$i}})/$nBinSize }
         elsif ($sType eq 'nfr')    { $nBinMeasure = &percentile($hBins{$i},10)    }
         else{ die "Error: unknown measurement type in get_measurement_bins()\n";  }
         push @asResult, $nBinMeasure;
      }
      else{
         if ( ($sType eq 'sum') or ($sType eq 'cov') ){
            push @asResult, 0;
         }
         else{
            push @asResult, 'NA';
         }
      }

   }
   return @asResult;
}


# mean
#
# Calculate mean of numeric array
sub mean {
    my(@data)=@_;
    my $sum;
    foreach(@data) {
	$sum+=$_;
    }
    return($sum/@data);
}


# median
#
# Calculate median of numeric array
sub median {
    my @data   = sort { $a <=> $b} @_;
    my $nCount = scalar(@data);
    if ($nCount % 2) {
	return($data[int($nCount/2)]);
    } 
    else {
	return(($data[$nCount/2] + $data[$nCount/2 - 1])/2);
    }
}


# sum
#
# Calculate sum of numeric array
sub sum {
    my(@data)=@_;
    my $sum;
    foreach(@data) {
	$sum+=$_;
    }
    return($sum);
}


# percentile
#
# Calculate percentile of numeric array
sub percentile {
  my $rArray      = shift;
  my $nPercentile = shift || 0;
  ##Since we're returning a single value there's no real need
  ##to cache this.

  ##If the requested percentile is less than the "percentile bin
  ##size" then return undef.  Check description of RFC 2330 in the
  ##POD below.
  my $nCount = scalar (@$rArray);
  return undef if $nPercentile < 100 / $nCount;

  my @anData = sort {$a<=>$b} @$rArray;
  my $num    = $nCount*$nPercentile/100;
  my $index  = &POSIX::ceil($num) - 1;
  return wantarray
    ?  ($anData[ $index ], $index)
    :   $anData[ $index ];
}


# print_relative
#
# Print relative coordinates instead of a big data matrix
sub print_relative {
   my ($sJoinOutput, $rFeaturesHash, $sMappingType, $nFlankSize) = @_;
   
   open JOIN, $sJoinOutput or die "Error: can't open join output file: $!\n";
   while (<JOIN>){
      next if (/^\s*$/);
      next if (/^ *#/);
      s/[\n\r]+$//;
      my ($sFt, $sJoinType, $nPos, $nVal) = split /;/;
      $nPos = $nPos - $nFlankSize if ($sJoinType eq 'u');
      $nPos = $nPos + $rFeaturesHash->{$sFt}{parent_spliced_size} if ( ($sJoinType eq 'd') and ($sMappingType eq 'feature') );
      print join("\t", $sFt, $rFeaturesHash->{$sFt}{strand},$rFeaturesHash->{$sFt}{type}, $nPos, $nVal), "\n";
   }
   close JOIN;
}


# INTERRUPT
#
# Interrupt routine, make sure we exit gracefully for tmp file cleanup
sub INTERRUPT{
   exit(1); # This will call END
}
