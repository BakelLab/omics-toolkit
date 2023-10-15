#!/usr/bin/env perl

# 08.01.2014 15:09:37 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use File::Temp qw(tempfile tempdir);

# GLOBALS
$ENV{TMPDIR}         ||= "/sc/hydra/scratch/$ENV{USER}/";   # location for tmp file storage
$ENV{BEDTLS}         ||= "bedtools";                         # location of bedtools binary
$SIG{'INT'}=$SIG{'HUP'}=$SIG{'ABRT'}=$SIG{'QUIT'}=$SIG{'TRAP'}=$SIG{'STOP'}=\&INTERRUPT;
my %hMultiCopyTags = ('tag'=>'','Dbxref'=>'','gbkey'=>'','Name'=>'','product'=>'','Note'=>'','tss_id'=>'', 'ont'=>'', 'db_xref'=>'', 'gene_synonym'=>'');

# GET PARAMETERS
my $sHelp        = 0;
my $nRegionSize  = 5000;
my $flThreePrime = 0;
my $flFivePrime  = 0;
my $nMinlength   = 500;
my $sRemoveOvlp  = "full";
my $flSameStrand = 0;
my $nBufferSize  = 0;
GetOptions("help!"        => \$sHelp,
           "regionsize:s" => \$nRegionSize,
           "bufferize:s"  => \$nBufferSize,
           "minlength:s"  => \$nMinlength,
           "3!"           => \$flThreePrime,
           "5!"           => \$flFivePrime,
           "s!"           => \$flSameStrand,
           "overlap:s"    => \$sRemoveOvlp);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [-3 -5] -r <length> <gtf-file>
   
   Extract coordinates of gene-flanking regions. For every gene
   the most distal coordinate of any annotated transcript is taken
   as the start coordinate and the end coordinate is extended 
   outward by the defined window size.
   
   Arguments:
    -r --regionsize <integer>
      Length of the feature-flanking region in bp
      Default: $nRegionSize
    -b --buffersize <integer>
      Optional length of a buffer region directly adjacent the 3' or 5' 
      end to exclude from the flanking region. Default: $nBufferSize
    -3
      Write 3' flanking regions
    -5
      Write 5' flanking regions
    -o --overlap <string>
      Optionally trim overlaps with other features. Options include
      'full' (remove any feature overlap), 'exon' (remove overlap with exons only),
      or 'none' (do not remove overlaps);
      Default: full
    -s
      Require same strandednes for overlap removal
    -m --minlength <integer>
      Minimum length of flanking region after trimming overlaps and chromosome ends.
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check arguments
$sRemoveOvlp = lc($sRemoveOvlp);
die "Error: Flanking region size must be an integer greater than 1.\n" unless ($nRegionSize =~ /^\d+$/);
die "Error: Buffer region size must be an integer.\n" unless ($nBufferSize =~ /^\d+$/);
die "Error: Buffer region size must be smaller than flanking region size" unless ($nBufferSize < $nRegionSize);
die "Error: Specify either -3 or -5 to extract 3' or 5' flanking regions.\n" unless($flThreePrime or $flFivePrime);
die "Error: Arguments -3 and -5 cannot be specified at the same time.\n" if($flThreePrime and $flFivePrime);
die "Error: Select either 'full', 'exon', or 'none' for the overlap argument.\n" unless($sRemoveOvlp =~ /^(full|exon|none)$/);

# Collect most distal start and end coordinates for each gene
my %hCoord;
open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]': $!\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^\s*#/);
   s/[\n\r]+$//;
   my ($sSeqName, $sSource, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sGroup) = split /\t/, $_, -1;
   ($nStart, $nEnd) = ($nEnd, $nStart) if ($nEnd < $nStart);
   my $rhAnnots = gtf_annots_to_hash($sGroup, \%hMultiCopyTags);
   if (exists $rhAnnots->{gene_id}){
      my $sGeneID = $rhAnnots->{gene_id};
      if ( exists $hCoord{$sGeneID} ){
         die "Error: gene '$sGeneID' has annotations on multiple strands"     if ($hCoord{$sGeneID}{strand} ne $sStrand);
         die "Error: gene '$sGeneID' has annotations on multiple chromosomes" if ($hCoord{$sGeneID}{chr}    ne $sSeqName);
         $hCoord{$sGeneID}{start}  = $nStart if ($nStart < $hCoord{$sGeneID}{start});
         $hCoord{$sGeneID}{end}    = $nEnd   if ($nEnd   > $hCoord{$sGeneID}{end});
      }
      else{
         $hCoord{$sGeneID}{start}  = $nStart;
         $hCoord{$sGeneID}{end}    = $nEnd;
         $hCoord{$sGeneID}{strand} = $sStrand;
         $hCoord{$sGeneID}{chr}    = $sSeqName;
         $hCoord{$sGeneID}{source} = $sSource;
      }
   }
}
close IN;

# Write flanking regions to a temporary file
my ($fhTmpFlank, $sTmpFlank) = tempfile('flanking-regions-unfiltered-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>0);
foreach my $sGeneID (keys %hCoord){
   my ($nOutStart, $nOutEnd) = (0,0);
   if ( ($hCoord{$sGeneID}{strand} eq "+") or ($hCoord{$sGeneID}{strand} eq "1") ){
      if ($flThreePrime){
         $nOutStart = $hCoord{$sGeneID}{end} + 1 + $nBufferSize;
         $nOutEnd   = $hCoord{$sGeneID}{end} + $nRegionSize + 1;
      }
      else{
         $nOutStart = $hCoord{$sGeneID}{start} - $nRegionSize - 1;
         $nOutEnd   = $hCoord{$sGeneID}{start} - 1 - $nBufferSize;
      }
   }
   else{
      if ($flThreePrime){
         $nOutStart = $hCoord{$sGeneID}{start} - $nRegionSize - 1;
         $nOutEnd   = $hCoord{$sGeneID}{start} - 1 - $nBufferSize;
      }
      else{
         $nOutStart = $hCoord{$sGeneID}{end} + 1 + $nBufferSize;
         $nOutEnd   = $hCoord{$sGeneID}{end} + $nRegionSize + 1;
      }
   }
   $nOutStart = 1 if ($nOutStart < 1);
   $nOutEnd   = 1 if ($nOutEnd < 1);
   if ( ($nOutEnd - $nOutStart) >= $nMinlength ){
      print $fhTmpFlank join("\t", $hCoord{$sGeneID}{chr}, $hCoord{$sGeneID}{source}, "exon", $nOutStart, $nOutEnd, ".", $hCoord{$sGeneID}{strand}, ".", "gene_id \"$sGeneID\";"), "\n";
   }
}
close($fhTmpFlank);

# Optionally remove overlap with other features
my $sStrandedness = $flSameStrand ? "-s" : "";
if ($sRemoveOvlp eq "full"){
   # Make a temporary feature file for bedtools
   my ($fhTmpFeature, $sTmpFeature) = tempfile('feature-regions-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>0);
   foreach my $sGeneID (keys %hCoord){
      print $fhTmpFeature join("\t", $hCoord{$sGeneID}{chr}, $hCoord{$sGeneID}{source}, "exon", $hCoord{$sGeneID}{start}, $hCoord{$sGeneID}{end}, ".", $hCoord{$sGeneID}{strand}, ".", "gene_id \"$sGeneID\";"), "\n";
   }
   close($fhTmpFlank);
   
   # Run bedtools subtract and output subtraction results
   open BEDTOOLS, "$ENV{BEDTLS} subtract -nonamecheck $sStrandedness -a $sTmpFlank -b $sTmpFeature |" or die "Error: can't run bedtools: $!\n";
   while (<BEDTOOLS>){
      next if (/^\s*$/);
      next if (/^ *#/);
      my ($sSeqName, $sSource, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sGroup) = split /\t/, $_, -1;
      print $_ if ( ($nEnd - $nStart) >= $nMinlength);
   }
   close BEDTOOLS;
}
elsif ($sRemoveOvlp eq "exon"){
   # Make a temporary feature file for bedtools
   my ($fhTmpFeature, $sTmpFeature) = tempfile('feature-regions-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>0);
   open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]': $!\n";
   while (<IN>){
      next if (/^\s*$/);
      next if (/^ *#/);
      my ($sSeqName, $sSource, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sGroup) = split /\t/, $_, -1;
      print $_ if ( $sFeature eq "exon");
   }
   close IN;
   close($fhTmpFlank);
   
   # Run bedtools bedtools subtract and output subtraction results
   open BEDTOOLS, "$ENV{BEDTLS} subtract -nonamecheck $sStrandedness -a $sTmpFlank -b $sTmpFeature |" or die "Error: can't run bedtools: $!\n";
   while (<BEDTOOLS>){
      next if (/^\s*$/);
      next if (/^ *#/);
      my ($sSeqName, $sSource, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $nFrame, $sGroup) = split /\t/, $_, -1;
      print $_ if ( ($nEnd - $nStart) >= $nMinlength);
   }
   close BEDTOOLS;
}
else{
   open OUT, $sTmpFlank or die "Error: can't open temporary file '$sTmpFlank':$!\n";
   while (<OUT>){ print }
   close OUT;
}


#################
## SUBROUTINES ##
#################


# gtf_annots_to_hash
#
# Parse gtf annotations and return key-value pairs
sub gtf_annots_to_hash {
   my ($sAnnots, $rMultiCopy) = @_;
   my %hReturn;
   
   my @asPairs = split /" *; */, $sAnnots;
   foreach my $sPair (@asPairs){
      my ($sKey, $sVal) = split / "/, $sPair;
      if ( exists $rMultiCopy->{$sKey} ){
         push @{$hReturn{$sKey}}, $sVal;
      }
      else{
         die "Error: Duplicate key entry '$sKey' found with values '$hReturn{$sKey}' and '$sVal'\n" if (exists $hReturn{$sKey});
         $hReturn{$sKey} = $sVal;
      }
   }
   return \%hReturn;
}

# INTERRUPT
#
# Interrupt routine, make sure we exit gracefully for tmp file cleanup
sub INTERRUPT{
   exit(1); # This will call END
}
