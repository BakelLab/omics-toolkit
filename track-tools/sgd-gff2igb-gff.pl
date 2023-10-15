#!/usr/bin/env perl

# Convert the SGD gff file into a gff file that can be used by IGB

# MODULES
use strict;
use Getopt::Long;

# GLOBALS
my %hGenes        = map { $_ => 0 } qw(CDS);
my %hNonCoding    = map { $_ => 0 } qw(tRNA snoRNA rRNA ncRNA snRNA);
my %hTransFacs    = map { $_ => 0 } qw(TF_binding_site);
my %hOther        = map { $_ => 0 } qw(long_terminal_repeat ARS transposable_element_gene LTR_retrotransposon 
                                       pseudogene telomere centromere internal_transcribed_spacer_region 
                                       external_transcribed_spacer_region);
my %hChrIDs       = (chrI   => "chr1",  chrII   => "chr2",  chrIII  => "chr3",  chrIV  => "chr4",  chrV  => "chr5",
                     chrVI  => "chr6",  chrVII  => "chr7",  chrVIII => "chr8",  chrIX  => "chr9",  chrX  => "chr10",
                     chrXI  => "chr11", chrXII  => "chr12", chrXIII => "chr13", chrXIV => "chr14", chrXV => "chr15",
                     chrXVI => "chr16", chrMito => "chr17", chrMT   => "chr17");
my %hTrackColors  = (sgd_orfs => "0,140,0", sgd_noncoding => "96,194,194", sgd_other => "190,190,190", sgd_transfac => "190,134,135");


# GET ARGUMENTS
my $sHelp     = 0;
my $sInput    = '';
GetOptions("help!"        => \$sHelp,
           "input:s"      => \$sInput);



# PRINT HELP
$sHelp = 1 unless($sInput);
if ($sHelp) {
    die <<HELP

    sgd-gff2igb-gff.pl -i input

    arguments:
    -i Input file
        Gff file downloaded from the SDG website (saccharomyces_cerevisiae.gff
        and scerevisiae_regulatory.gff)
    -h Help

HELP
}

###########
## START ##
###########


# Open file and read features into one of four arrays (Genes, Non-coding, Other, Transfacs)
my %hResult;
my %hSkippedFeatures;
open INPUT, $sInput or die "Can't open sgd gff file: $!\n";
while (<INPUT>){
   next if (/^\s*$/);
   s/[\n\r]//g;
   my ($chr,$source,$feature,$start,$end,$score,$strand,$frame,$group) = split /\t/;
   my $id  = &get_identifier_from_gff3_group($group);
   $score  = 1000;
   $chr    = $hChrIDs{$chr} if (exists($hChrIDs{$chr}));
   
   # Divide the different features in different categories
   if    ( exists($hGenes{$feature}) ){
      my $sOrfStatus = &get_orf_status($group);
      $score = 500 if (lc($sOrfStatus) eq 'uncharacterized');
      $score = 4   if (lc($sOrfStatus) eq 'dubious');
      
      # Add ORF to output, only if the status is defined (otherwise it's a pseudogene, or transposable element).
      push @{$hResult{sgd_orfs}}, [$chr,$source,$feature,$start,$end,$score,$strand,$frame,$id] if ($sOrfStatus);
   }
   elsif ( exists($hNonCoding{$feature}) ){
      push @{$hResult{sgd_noncoding}}, [$chr,$source,$feature,$start,$end,$score,$strand,$frame,$id];
   }
   elsif ( exists($hOther{$feature}) ) {
      push @{$hResult{sgd_other}}, [$chr,$source,$feature,$start,$end,$score,$strand,$frame,$id];
   }
   elsif ( exists($hTransFacs{$feature}) ){
      push @{$hResult{sgd_transfac}}, [$chr,$source,$feature,$start,$end,$score,$strand,$frame,$id];
   }
   else{
      $hSkippedFeatures{$feature}++;
   }
      
}


# Print a track file for each category
foreach my $sCategory (keys(%hResult)){
   open OUTPUT, ">$sCategory.gff" or die "Can't write to $sCategory.gff: $!\n";
   print OUTPUT "##gff-version\t2\n";
   print OUTPUT "track name=\"$sCategory\" useScore=1 color=$hTrackColors{$sCategory} url=\"http://db.yeastgenome.org/cgi-bin/locus.pl?locus=\$\$\"\n";
   foreach my $rKey (@{$hResult{$sCategory}}){
      print OUTPUT join("\t", @$rKey), "\n";
   }
   close OUTPUT;
}


#################
## SUBROUTINES ##
#################



# get_identifier_from_gff3_group(group_string)
#
# Extracts a single identifier from the structured gff3-formatted group field
# Extraction order is: Gene, Name, ID
sub get_identifier_from_gff3_group{
   my $sGffGroup = shift @_;
   my $sReturn = '';
   
   # Split group into hash
   my %hGffGroup;
   my @asPairs = split /;/, $sGffGroup;
   foreach my $sPair (@asPairs){
      my ($key, $val) = split /=/, $sPair;
      $hGffGroup{lc($key)} = $val;
   }
   
   # Assign key in order of preference; gene, name, id
   my $sID;
   if     ( exists($hGffGroup{gene}) ) { $sID = $hGffGroup{gene} }
   elsif  ( exists($hGffGroup{name}) ) { $sID = $hGffGroup{name} }
   elsif  ( exists($hGffGroup{id}) )   { $sID = $hGffGroup{id}   }
   else                                { $sID = 'Unknown'        }

   # Some custom post-processing steps to pretty up the IDs
   $sID = &parse_identifier($sID);
   
   # format a gff2 group
   my @asReturn;
   push @asReturn, "gene_id \"$sID\"";
   push @asReturn, "transcript_id \"$sID\"";
   if (exists($hGffGroup{note})){
      my $sNote = &parse_identifier($hGffGroup{note});
      push @asReturn, "note \"$sNote\"";
   }
   if (exists($hGffGroup{orf_classification})){
      my $sClassification = &parse_identifier($hGffGroup{orf_classification});
      push @asReturn, "orf_classification \"$sClassification\"";
   }
   
   # Return ID
   $sReturn = join(';', @asReturn);
   return($sReturn);
   
}


# parse_identifier
#
# Remove formatting codes from string
sub parse_identifier{
   my $sID = shift @_;
   $sID =~ s/%20/ /g;
   $sID =~ s/%24/\$/g;
   $sID =~ s/%26/\&/g;
   $sID =~ s/%2B/\+/g;
   $sID =~ s/%2C/\,/g;
   $sID =~ s/%2F/\//g;
   $sID =~ s/%3A/\:/g;
   $sID =~ s/%3B/\./g;
   $sID =~ s/%3D/\=/g;
   $sID =~ s/%3F/\?/g;
   $sID =~ s/%40/\@/g;
   $sID =~ s/%25/\%/g;
   $sID =~ s/%5B/\[/g;
   $sID =~ s/%5D/\]/g;
   $sID =~ s/%3E/\>/g;
   $sID =~ s/%22/\'/g;
   return ($sID);
}


# get_orf_status
#
# Gets the status of an ORF: verified, dubious or uncharacterized
# Returns an empty value if the status is not defined.
sub get_orf_status{
   my $sGffGroup = shift @_;
   my $sReturn = '';
   
   # Split group into hash
   my %hGffGroup;
   my @asPairs = split /;/, $sGffGroup;
   foreach my $sPair (@asPairs){
      my ($key, $val) = split /=/, $sPair;
      $hGffGroup{lc($key)} = $val;
   }
   
   # Get status
   $sReturn = lc($hGffGroup{orf_classification}) if exists($hGffGroup{orf_classification});
   
   return ($sReturn);
}