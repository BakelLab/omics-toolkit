#!/usr/bin/env perl

# Takes a bpmap file and extracts the specified probe groups into a
# tab-delimited file that can directly be used by the 
# map-array-probes_run-blat script.

# MODULES
use strict;
use Getopt::Long;
use File::Basename;

# ENVIRONMENT
$ENV{'BPMAP2TPMAP_BINARY'} ||= 'bpmap2tpmap';
$ENV{'QSUB_BINARY'}        ||= 'qsub';

# GET ARGUMENTS
my $sHelp     = 0;
my $sInput    = '';
my $sSeqGroup = '';
GetOptions("help!"      => \$sHelp,
	   "input:s"    => \$sInput,
	   "seqgroup:s" => \$sSeqGroup);

# PRINT HELP
$sHelp = 1 unless($sInput and $sSeqGroup);
if ($sHelp) {
    die <<HELP
    
    map-array-probes_extract-bpmap-probes -i <bpmap file> -s <sequence group name>
    
    Takes a bpmap file and extracts the specified probe groups into a
    tab-delimited file that can directly be used by the 
    map-array-probes_run-blat script. A unique identifier is created for 
    each probe, based on the PM position on the array. The information in 
    this identifier is later used when parsing the BLAT output and creating 
    the remapped bpmap file.
    
    arguments:
    -i Bpmap file
        Tab delimited list of probe ID and sequence
         probeID <tab> sequence
    -s Sequence group name
        The name of the sequence groups to extract from the bpmap file
        Only these probes will be remapped to the new genome version.
    -h Help
    
HELP
}


######################################
#               START                #
######################################

my %hCheck;
open BPMAP, "$ENV{'BPMAP2TPMAP_BINARY'} -bpmap $sInput -group $sSeqGroup |" or die "Can't extract bpmap file: $!\n";
while (<BPMAP>){
   if (/^#probeset_type/){
      s/\n|\r//g;
      my ($sKey, $sValue) = split /\t/;
      die "Incompatible probe set type '$sValue' found. Can only process tiling probes!\n" unless ($sValue eq 'tiling');
   }
   next if (/^\s*$/);
   next if (/^\s*#/);
   s/\n|\r//g;
   my (@asLine) = split /\t/;
   my $sID      = join('|', $sSeqGroup, @asLine[(4,5)]);
   print join("\t", $sID, $asLine[0]), "\n";
   $hCheck{$sID}++;
}

# Check to make sure that there are no duplicate keys
my $flCheck = 0;
foreach my $sKey (keys(%hCheck)){
   if ($hCheck{$sKey}>1){
      $flCheck = 1;
      last;
   }
}
if ($flCheck){
   print STDERR "WARNING: Duplicate probe identifiers found, file should not be used for remapping!!\n";
}
else{
   print STDERR "No problems found, proceed with BLAT analysis\n";
}
