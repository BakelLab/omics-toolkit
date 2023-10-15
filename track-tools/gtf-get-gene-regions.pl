#!/usr/bin/env perl

# 08.01.2014 15:09:37 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use File::Temp qw(tempfile tempdir);

# GLOBALS
$ENV{TMPDIR}         ||= "/sc/orga/scratch/vanbah01/tmp/";   # location for tmp file storage
$ENV{BEDTLS}         ||= "bedtools";                         # location of bedtools binary
$SIG{'INT'}=$SIG{'HUP'}=$SIG{'ABRT'}=$SIG{'QUIT'}=$SIG{'TRAP'}=$SIG{'STOP'}=\&INTERRUPT;
my %hMultiCopyTags = ('tag'=>'','Dbxref'=>'','gbkey'=>'','Name'=>'','product'=>'','Note'=>'','tss_id'=>'', 'ont'=>'', 'db_xref'=>'', 'gene_synonym'=>'');

# GET PARAMETERS
my $sHelp        = 0;
my $n3PrimeTrim  = 0;
my $n5PrimeTrim  = 0;
my $nMinlength   = 500;
GetOptions("help!"        => \$sHelp,
           "minlength:s"  => \$nMinlength,
           "3:s"           => \$n3PrimeTrim,
           "5:s"           => \$n5PrimeTrim);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [-3 -5] -r <length> <gtf-file>
   
   Extract coordinates of gene regions. For every gene
   the most distal coordinate of any annotated transcript is taken
   as the start coordinate and the end coordinate. Starts and
   ends can be trimmed by an optional defined length. The score
   field will contain the maximum exon count for gene transcripts
   and can be used to differentiate between e.g. mono-exon and
   multi-exon genes.
   
   Arguments:
    -3 <integer>
      Length to trim back from the 3' end. Default: $n3PrimeTrim
    -5 <integer>
      Length to trim back from the 5' end. Default: $n5PrimeTrim
    -m --minlength <integer>
      Minimum length of gene regions after trimming ends. Default: $nMinlength
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check arguments
die "Error: -3 argument must be an integer.\n" unless($n3PrimeTrim =~ /^\d+$/);
die "Error: -5 argument must be an integer.\n" unless($n5PrimeTrim =~ /^\d+$/);

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
      my $sTranscriptID = $rhAnnots->{transcript_id};
      if ( exists $hCoord{$sGeneID} ){
         die "Error: gene '$sGeneID' has annotations on multiple strands"     if ($hCoord{$sGeneID}{strand} ne $sStrand);
         die "Error: gene '$sGeneID' has annotations on multiple chromosomes" if ($hCoord{$sGeneID}{chr}    ne $sSeqName);
         $hCoord{$sGeneID}{start}  = $nStart if ($nStart < $hCoord{$sGeneID}{start});
         $hCoord{$sGeneID}{end}    = $nEnd   if ($nEnd   > $hCoord{$sGeneID}{end});
         $hCoord{$sGeneID}{exons}{$sTranscriptID}++ if (lc($sFeature) eq 'exon');
      }
      else{
         $hCoord{$sGeneID}{start}  = $nStart;
         $hCoord{$sGeneID}{end}    = $nEnd;
         $hCoord{$sGeneID}{strand} = $sStrand;
         $hCoord{$sGeneID}{chr}    = $sSeqName;
         $hCoord{$sGeneID}{source} = $sSource;
         $hCoord{$sGeneID}{exons}{$sTranscriptID}++ if (lc($sFeature) eq 'exon');
      }
   }
}
close IN;

# Write gene regions
foreach my $sGeneID (keys %hCoord){
   my ($nOutStart, $nOutEnd) = (0,0);
   
   # Consider flanking
   if ( ($hCoord{$sGeneID}{strand} eq "+") or ($hCoord{$sGeneID}{strand} eq "1") ){
      $nOutStart = $hCoord{$sGeneID}{start} + $n5PrimeTrim;
      $nOutEnd   = $hCoord{$sGeneID}{end} - $n3PrimeTrim;
   }
   else{
      $nOutStart = $hCoord{$sGeneID}{start} + $n3PrimeTrim;
      $nOutEnd   = $hCoord{$sGeneID}{end} - $n5PrimeTrim;
   }
   
   # Get exon count
   my $nMaxExonCount = 0;
   foreach my $sTranscriptID ( keys %{$hCoord{$sGeneID}{exons}} ){
      $nMaxExonCount = $hCoord{$sGeneID}{exons}{$sTranscriptID} if ( $hCoord{$sGeneID}{exons}{$sTranscriptID} > $nMaxExonCount) ;
   }
   
   if ( ($nOutEnd - $nOutStart) >= $nMinlength ){
      print join("\t", $hCoord{$sGeneID}{chr}, $hCoord{$sGeneID}{source}, "exon", $nOutStart, $nOutEnd, $nMaxExonCount, $hCoord{$sGeneID}{strand}, ".", "gene_id \"$sGeneID\";"), "\n";
   }
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
