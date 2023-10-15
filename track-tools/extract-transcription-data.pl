#!/usr/bin/env perl

# extract-transcription-data.pl
# Extracts the data from the David paper on whole genome transcription in yeast from the paper

# MODULES

use strict;
use Getopt::Long;
use FileHandle;
use Bio::Seq;
use Bio::SeqIO;
use Cwd;
require "/home/gen/harm/script/igb/igb-conf.pm";

# GET ARGUMENTS
my $sHelp     = 0;
my $sInput    = '';
my $sChrDir   = '';
my $nOffset   = 25;
GetOptions("help!"    => \$sHelp,
	   "input:s"  => \$sInput,
	   "dir:s"    => \$sChrDir,
	   "offset:s" => \$nOffset
	   );

# PRINT HELP
$sHelp = 1 unless($sInput and $sChrDir);
if ($sHelp) {
    die <<HELP
    extract-transcription-data.pl
    Maps whole-genome transcription data
    
    arguments:
    -i Input file
        Input formatted according to the original paper.
    -d Chromosome directory
        Directory containing original chromosome sequence files with names
        formatted as 'chrNN.fsa'. There are no leading zeros for the chromosome
        name.
    -o Offset
        Amount of sequence that is appended before and after the original sequence
        motif. This extended sequence containing the motif will be used for blasting
        and subsequent mapping of the motif to a new igb genome.
        Default: 25
    -h Help
HELP
}


######################################
#               START                #
######################################

unless (-e $sChrDir) {die "Can't read from '$sChrDir': directory does not exist"};
my $rhBindingSites = &get_gff_hash($sInput);
foreach my $sChromosome (keys %{$rhBindingSites}){
    my $sFilePath = join('/', $sChrDir, "$sChromosome.fsa");
    die "Can't read from '$sFilePath': file does not exist" unless (-e $sFilePath);
    my $oFasta    = Bio::SeqIO->new(-file => $sFilePath , '-format' => 'fasta'); 
    my $oSeq      = $oFasta->next_seq();
    
    foreach my $rPosition (@{$rhBindingSites->{$sChromosome}}){
	my ($sStrand, $nPosition, $nPolA, $nTotal) = @$rPosition;
	my $nStart       = $nPosition - $nOffset;
	$nStart       = 1 if ($nStart < 1);
	my $nStop        = $nPosition + $nOffset;
	$nStop        = $oSeq->length() if ($nStop > $oSeq->length());
	my $nLength  = $nStop - $nStart + 1;
	my $sSeq = $oSeq->subseq($nStart,$nStop);
	print ">$sChromosome|$nPosition|$nPolA|$sStrand|$nLength\n$sSeq\n";
    }
}


######################################
#            SUBROUTINES             #
######################################

# Reads the input file and returns a hash with values for each chromosome

sub get_gff_hash {
    my $sFile = shift @_;
    
    my %hBindingSites;
    my $fhInput = new FileHandle;
    $fhInput->open($sFile) or die "Can't open input file: $!\n";
    while (<$fhInput>){
	next if (/^\s*$/);
	next if (/NA/);
	s/[\n\r]//g;
	my ($sChr, @asLine) = split /\t/;
	
	# some checks
	unless ($sChr =~ 'chr'){
	    print STDERR "Warning: chromosome names not formatted as 'chrNN', please remap to allow loading in IGB.\n"; }
	unless (scalar(@asLine) == 4){
	    die "Error: the number of columns in the file is not consistent with the expected file format\n"; }
	
	# put the stuff in the hash
	push @{$hBindingSites{$sChr}}, [@asLine];
    }
    return (\%hBindingSites);
}


