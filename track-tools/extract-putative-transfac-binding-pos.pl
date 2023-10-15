#!/usr/bin/env perl

# extract-putative-transfac-binding-pos.pl
# Extracts putative transcription factor binding sites from the rick young data

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
my $nOffset   = 20;
GetOptions("help!"    => \$sHelp,
	   "input:s"  => \$sInput,
	   "dir:s"    => \$sChrDir,
	   "offset:s" => \$nOffset
	   );

# PRINT HELP
$sHelp = 1 unless($sInput and $sChrDir);
if ($sHelp) {
    die <<HELP
    extract-putative-transfac-binding-pos.pl
    To map Rick young data
    
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
	my ($sSite, $nStart, $nStop, $sStrand) = @$rPosition;
	$nStart       = $nStart - $nOffset;
	$nStart       = 1 if ($nStart < 1);
	$nStop        = $nStop + $nOffset;
	$nStop        = $oSeq->length() if ($nStop > $oSeq->length());
	my $nLength  = $nStop - $nStart + 1;
	my $sSeq = $oSeq->subseq($nStart,$nStop);
	print ">$sSite|$sChromosome|$nStart|$nStop|$nOffset|$nLength|$sStrand\n$sSeq\n";
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
	s/[\n\r]//g;
	my ($sChr, @asLine) = split /\t/;
	
	# some checks
	unless ($sChr =~ 'chr'){
	    print STDERR "Warning: chromosome names not formatted as 'chrNN', please remap to allow loading in IGB.\n"; }
	unless (scalar(@asLine) == 8){
	    die "Error: the number of columns in the file is not consistent with the expected file format\n"; }
	
	# put the stuff in the hash
	my ($sSite, $sStart, $sStop, $sStrand) = ($asLine[7],$asLine[2], $asLine[3], $asLine[5]);
	$sSite =~ s/Site //g;
	$sSite =~ s/;//g;
	push @{$hBindingSites{$sChr}}, [$sSite, $sStart, $sStop, $sStrand];
    }
    return (\%hBindingSites);
}


