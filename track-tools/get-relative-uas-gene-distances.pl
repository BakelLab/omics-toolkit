#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Bio::SeqFeature::Collection;
use Bio::SeqFeature::Generic;
use Bio::Graphics;


# GET ARGUMENTS
my $sHelp       = 0;
my $sTransfacs  = '';
my $sGenes      = '';
my $nPromSize   = 800;
GetOptions("help!"       => \$sHelp,
	   "transfacs:s" => \$sTransfacs,
	   "genes:s"     => \$sGenes,
	   "prom:s"      => \$nPromSize);

# PRINT HELP
$sHelp = 1 unless($sTransfacs and $sGenes and $nPromSize);
if ($sHelp) {
    die <<HELP
    get-relative-uas-gene-distances.pl -t <transfac.gff> -g <genes.gff>
    Returns a list of genes with a variety of statistics on the positions of UAS
    
    arguments:
    -t Transcription factor list
        List of transcription factors in GFF format
    -g Gene list
        List of genes in peakcluster format (IGB directory)
    -p Promotor size
        Default: 800
    -h Help
HELP
}


######################################
#               START                #
######################################

# Read the transcription factor list and put in hash
my %hTF;
open INPUT, $sTransfacs or die "Can't open transcription factor file '$sTransfacs': $!";
while (<INPUT>){
    next if (/^\s*$/);
    s/[\n|\r]//g;
    my ($sChr, $sSource, $sType, $nStart, $nStop, $nScore, $nStrand, $nBla, $sName) = split /\t/;
    my $oBioFeature   = Bio::Graphics::Feature->new(-display_name => $sName,
						    -start        => $nStart,
						    -end          => $nStop);
    push @{$hTF{$sChr}}, $oBioFeature;
}
close INPUT;


# Convert hash to annotation object per chromosome
foreach my $sChr (sort(keys(%hTF))){
    my @aoFeatures  = @{$hTF{$sChr}};
    my $oCollection = new Bio::SeqFeature::Collection();
    my $nAdded      = $oCollection->add_features(\@aoFeatures);
    $hTF{$sChr}     = $oCollection;
    print STDERR "added $nAdded transcription factor binding sites for chromosome $sChr\n";
}


# Hash the genelist to make sure that multi-exon genes are grouped together
my %hGenes;
open GENES, $sGenes or die "Can't open gene list '$sGenes': $!";
while (<GENES>){
    next if (/^\s*$/);
    s/[\n|\r]//g;
    my ($sID, $sType, $sChr, $nStart, $nStop, $sStrand) = split /\t/;
    if (exists($hGenes{$sID})){
	my ($nTmpStart, $nTmpStop, $sTmpStrand, $sTmpChr) = @{$hGenes{$sID}};
	$nTmpStart = $nStart if ($nStart < $nTmpStart);
	$nTmpStop  = $nStop  if ($nStop  > $nTmpStop);
	$hGenes{$sID} = [$nTmpStart, $nTmpStop, $sTmpStrand, $sTmpChr];
    }
    else{
	$hGenes{$sID} = [$nStart, $nStop, $sStrand, $sChr];
    }
} 
close GENES;


# Now get transcription factor positions for all genes
print "ID\tUAS-start\tUAS-stop\tUAS-mean\tUAS-median\n";
foreach my $sGene (sort(keys(%hGenes))){
    my ($nStart, $nStop, $sStrand, $sChr) = @{$hGenes{$sGene}};
    my ($nIntStart, $nIntStop);
    
    if ($sStrand eq '+'){
	$nIntStart = $nStart-$nPromSize;
	$nIntStart = 1 if ($nIntStart<1);
	$nIntStop  = $nStart;
    }
    else{
	$nIntStart = $nStop;
	$nIntStop  = $nStop+$nPromSize;
    }
    
    if (exists($hTF{$sChr})){
	my $oCollection = $hTF{$sChr};
	my @aoFeatures  = $oCollection->features_in_range(-start       => $nIntStart,
				             		  -end         => $nIntStop,
				             		  -strand      => 1,
				             		  -contain     => 0,
				             		  -strandmatch => 'ignore');
	if (@aoFeatures){
	    my @anTFdistances;
	    foreach my $oFeature (@aoFeatures){
		my $nPosition = int( ($oFeature->start()+$oFeature->end())/2 );
		my $nDistance;
		if ($sStrand eq '+') { $nDistance = $nIntStop  - $nPosition; }
		else                 { $nDistance = $nPosition - $nIntStart; }
		push @anTFdistances, $nDistance;
	    }
	    @anTFdistances = sort {$a<=>$b} @anTFdistances;
	
	    print join("\t", $sGene, $anTFdistances[0], $anTFdistances[$#anTFdistances], int(&mean(@anTFdistances)), int(&median(@anTFdistances))), "\n";
	}
    }
    else{
	print STDERR "Skipping gene $sGene: No transcription factor binding sites mapped to chromosome $sChr\n";
    }
}



#########################
#      SUBROUTINES      #
#########################


# Calculate mean of numeric array
sub mean {
    my(@data)=@_;
    my $sum;
    foreach(@data) {
	$sum+=$_;
    }
    return($sum/@data);
}

# Calculate median of numeric array
sub median {
    my(@data)=sort { $a <=> $b} @_;
    if (scalar(@data)%2) {
	return($data[@data/2]);
    } else {
	my($upper, $lower);
	$lower=$data[@data/2];
	$upper=$data[@data/2 - 1];
	return(mean($lower, $upper));
    }
}
