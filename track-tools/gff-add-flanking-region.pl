#!/usr/bin/env perl

# 15.08.2009 14:05:48 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use FileHandle;

# GET PARAMETERS
my $sHelp                = 0;
my $sInput               = '';
my $nFlankSize           = 500;
my $nFlankBufferSize     = 0;
my $sExcludeFlankOverlap = 'none';
my $flFlankOnly          = 0;
GetOptions("help!"                   => \$sHelp,
           "gff:s"                   => \$sInput,
           "flank_size:n"            => \$nFlankSize,
           "flank_buffer_size"       => \$nFlankBufferSize,
           "exclude_flank_overlap:s" => \$sExcludeFlankOverlap,
           "flanks_only!"            => \$flFlankOnly);

# PRINT HELP
$sHelp = 1 unless($sInput and $nFlankSize and $sExcludeFlankOverlap);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName
    
    -gff <string>
      The GFF file
    -flank_size <integer>
      The size of the upstream and downstream flanking regions.
      default: $nFlankSize
    -flank_buffer_size <integer>
      Additional buffer to increase distance to neighboring features
      default: $nFlankBufferSize
    -exclude_flank_overlap
      Prevent flanking regions from overlapping neighboring features.
       both      :   Prevent overlap with features on either strand
       sense     :   Prevent overlap with features on the same strand
       antisense :   Prevent overlap with features on the opposite strand
       none      :   Ignore overlapping features
      default: $sExcludeFlankOverlap
    -flanks_only
      Specify if you only want to output flanking regions (each on a separate line)
    -help
      This help message
   
HELP
}


###########
## START ##
###########

# Check arguments
$sExcludeFlankOverlap = lc($sExcludeFlankOverlap);
die "Error: Unknown flank overlap type '$sExcludeFlankOverlap'\n" unless($sExcludeFlankOverlap =~ /^(both|sense|antisense|none)$/);
die "Error: Flank size must be a discrete number\n"               unless($nFlankSize =~ /^\d+$/);
die "Error: Flank buffer size must be a discrete number\n"        unless($nFlankBufferSize =~ /^\d+$/);

# Read feature file
my %hFeatures;
read_gff_features(\%hFeatures, $sInput);

# Add start and end of closest flanking feature, if ExcludeOverlap is set
set_flank_positions(\%hFeatures, $nFlankSize, $nFlankBufferSize, $sExcludeFlankOverlap) unless ($sExcludeFlankOverlap eq 'none');

# Write output with flank size
foreach my $sFeature (keys %hFeatures){
   my $nFlankStart = $hFeatures{$sFeature}{'parent_start'} < $nFlankSize ? 0 : $hFeatures{$sFeature}{'parent_start'} - $nFlankSize;
   my $nFlankEnd   = $hFeatures{$sFeature}{'parent_end'} + $nFlankSize;
   if (exists $hFeatures{$sFeature}{'flank_start'}){
      $nFlankStart = $hFeatures{$sFeature}{'flank_start'} if ($hFeatures{$sFeature}{'flank_start'} > $nFlankStart);
   }
   if (exists $hFeatures{$sFeature}{'flank_end'}){
      $nFlankEnd   = $hFeatures{$sFeature}{'flank_end'}   if ($hFeatures{$sFeature}{'flank_end'} < $nFlankEnd);
   }
   
   # Write the output
   if ($flFlankOnly){
      print join("\t", $hFeatures{$sFeature}{'chr'}, 'perl', $hFeatures{$sFeature}{'type'}, $nFlankStart, $hFeatures{$sFeature}{'parent_start'}-1, 
                       '.', $hFeatures{$sFeature}{'strand'}, '.', $sFeature), "\n";
      print join("\t", $hFeatures{$sFeature}{'chr'}, 'perl', $hFeatures{$sFeature}{'type'}, $hFeatures{$sFeature}{'parent_end'}+1, $nFlankEnd, 
                       '.', $hFeatures{$sFeature}{'strand'}, '.', $sFeature), "\n";
   }
   else{
      print join("\t", $hFeatures{$sFeature}{'chr'}, 'perl', $hFeatures{$sFeature}{'type'}, $nFlankStart, $nFlankEnd, '.', 
                       $hFeatures{$sFeature}{'strand'}, '.', $sFeature), "\n";
   }
}


#################
## SUBROUTINES ##
#################



# read_gff_features($sFeatureFile)
#
# Reads the gff features into a hash
# returns hash with:  -> chromosome          : Chromosome
#                     -> strand              : Strand (+ or - or .)
#                     -> parent_start        : Start of first subfeature
#                     -> parent_end          : End of last subfeature
sub read_gff_features{
   my ($rhFeatures, $sFeaturesFile) = @_;
   
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
      
      # Make sure that we can extract an identifier from each line and do some other basic checks on gff parsing
      $sID = &extract_gff2_id($sID) if ($nGffVersion == 2);
      die "Error: line $. in '$sFeaturesFile' is missing an identifier (expecting GFF format)\n" unless $sID;
      die "Error: unknown strand '$sStrand' for feature '$sID', line $. in '$sFeaturesFile'\n" unless ($sStrand =~ /[+-.]/);
      die "Error: expected number for start position on line $. in '$sFeaturesFile'\n" unless ($nStart =~ /^\d+$/);
      die "Error: expected number for end position on line $. in '$sFeaturesFile'\n" unless ($nEnd =~ /^\d+$/);
      
      # Load the hash, note that we'll skip the array of features if flanking_only is set
      ($nStart, $nEnd) = ($nEnd, $nStart) if ($nStart > $nEnd); # Make sure that the start is always smaller than the end
      if (exists($rhFeatures->{$sID})) {
         die "Error: feature '$sID' maps to multiple chromosomes in the gff file\n"        unless ($rhFeatures->{$sID}{chr} eq $sChr);
         die "Error: feature '$sID' has mappings to both the forward and reverse strand\n" unless ($rhFeatures->{$sID}{strand} eq $sStrand);   
         $rhFeatures->{$sID}{parent_start} = $nStart if ($nStart < $rhFeatures->{$sID}{parent_start});
         $rhFeatures->{$sID}{parent_end}   = $nEnd if ($nEnd > $rhFeatures->{$sID}{parent_end});
      }
      else{
         $rhFeatures->{$sID}{chr}                 = $sChr;
         $rhFeatures->{$sID}{strand}              = $sStrand;
         $rhFeatures->{$sID}{type}                = $sType;
         $rhFeatures->{$sID}{parent_start}        = $nStart;
         $rhFeatures->{$sID}{parent_end}          = $nEnd;
      }
   }
   $fhGff->close();
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
            last if (  $aaPositions[$j][1] ne $sChr );                                            # only consider the same chromosome            
            next if ( ($aaPositions[$j][2] ne $sStr) and ($sExcludeFlankStrand eq 'sense') );     # only consider same strand, if specified
            next if ( ($aaPositions[$j][2] eq $sStr) and ($sExcludeFlankStrand eq 'antisense') ); # only consider opposite strand, if specified
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
            next if ( ($aaPositions[$j][2] ne $sStr) and ($sExcludeFlankStrand eq 'sense'));     # only consider same strand, if specified
            next if ( ($aaPositions[$j][2] eq $sStr) and ($sExcludeFlankStrand eq 'antisense')); # only consider opposite strand, if specified
            $nFlankEnd = $aaPositions[$j][3]-$nFlankBufferSize  if ($aaPositions[$j][3]-$nFlankBufferSize<$nFlankEnd);
         }
      }
      $nFlankEnd = $nEnd if ($nFlankEnd<$nEnd);                             # for overlapping features
      $rhFeatures->{$sFt}{flank_end} = $nFlankEnd;
   }  
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
