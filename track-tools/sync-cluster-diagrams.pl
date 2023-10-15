#!/usr/bin/env perl

# 25.01.2010 09:17:16 EST
# Harm van Bakel <hvbakel@gmail.com>

# GLOBALS
$ENV{TMPDIR}             ||= '/tmp';     # location for tmp file storage
$ENV{INTERSECT_ROW_IDS}  ||= 'sorted-intersect-by-ids';
$ENV{INTERSECT_COL_IDS}  ||= 'intersect-columns-by-ids';
$SIG{'INT'}=$SIG{'HUP'}=$SIG{'ABRT'}=$SIG{'QUIT'}=$SIG{'TRAP'}=$SIG{'STOP'}=\&INTERRUPT;

# MODULES
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Temp qw(tempfile tempdir);
use File::Copy;
 
# GET PARAMETERS
my $sHelp          = 0;
my $sSyncFile      = "";
my $sToFile        = "";
my $sOut       = "";
my $sArrayMappings = "";
my $sGeneMappings  = "";
my $sSyncOn      = "rows-and-columns";
GetOptions("help!"     => \$sHelp,
           "sync:s"    => \$sSyncFile,
           "to:s"      => \$sToFile,
           "out:s"     => \$sOut,
           "colmap:s"  => \$sArrayMappings,
           "rowmap:s"  => \$sGeneMappings,
           "on:s"      => \$sSyncOn);

# PRINT HELP
$sHelp = 1 unless($sSyncFile and $sToFile and $sOut);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [options] -sync <cdt file> -to <cdt file> -out <output file>
    
    -sync <string>
      The name of the .cdt file of to be synced
    -to <string>
      The name of the .cdt file to be synced to
    -on <string>
      The dimensions that the cluster diagrams should be synced on
      can be either 'rows' 'columns' or 'rows-and-columns'
      default: $sSyncOn
    -out <string>
      The basename of the output file
    -colmap <string>
      Optional file with column/array ID mappings between the 'sync' and 'to' file
      File format: id_to_file <tab> id_sync_file
    -rowmap <string>
      Optional file with row/gene ID mappings between the 'sync' and 'to' file
      File format: id_to_file <tab> id_sync_file
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Check Arguments
my %hSyncOn = ('rows'=>0,'columns'=>0,'rows-and-columns'=>0);
die "Error: unknown selection for the -on argument\n" unless(exists $hSyncOn{$sSyncOn});

# Get the line and column IDs for the relevant lines
my %hCdtPar;
$hCdtPar{sync} = get_cdt_params($sSyncFile, $sSyncOn);
$hCdtPar{to}   = get_cdt_params($sToFile, $sSyncOn);

# Start with syncing the gene IDs
my $sStagingFile = $sSyncFile;
if ( ($sSyncOn eq "rows") or ($sSyncOn eq "rows-and-columns") ){
   my $sOptions = "-ff $sStagingFile -fc $hCdtPar{sync}{gene_id_col} -fo $hCdtPar{sync}{gene_offset} ";
   $sOptions   .= "-if $sToFile      -ic $hCdtPar{to}{gene_id_col}   -io $hCdtPar{to}{gene_offset} ";
   $sOptions   .= "-mf $sGeneMappings" if ($sGeneMappings);
   
   my ($fhTmp, $sTmp) = tempfile('sync-clust-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1);
   open GENEFILT, "$ENV{INTERSECT_ROW_IDS} $sOptions |" or die "Error: gene ID filtering failed: $!\n";
   while (<GENEFILT>){
      if ($. >= $hCdtPar{sync}{gene_offset}){
         my @asLine = split /\t/, $_, -1;
         $asLine[0] = shift @{$hCdtPar{to}{genetree_ids}};
         $fhTmp->print(join("\t", @asLine));
      }
      else{
         $fhTmp->print($_);
      }
   }
   close GENEFILT;
   close $fhTmp;
   
   die "Error: could not sync files because '$sSyncFile' is missing geneIDs relative to '$sToFile'\n" if(scalar(@{$hCdtPar{to}{genetree_ids}}));
   $sStagingFile = $sTmp;
}

# Now do the same for the columns
if ( ($sSyncOn eq "columns") or ($sSyncOn eq "rows-and-columns") ){
   my $sOptions = "-ff $sStagingFile -fr $hCdtPar{sync}{array_id_row} -fo $hCdtPar{sync}{array_offset} ";
   $sOptions   .= "-if $sToFile      -ir $hCdtPar{to}{array_id_row}   -io $hCdtPar{to}{array_offset} -ro ";
   $sOptions   .= "-mf $sArrayMappings" if ($sArrayMappings);
   
   my ($fhTmp, $sTmp) = tempfile('sync-clust-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1);
   open ARRAYFILT, "$ENV{INTERSECT_COL_IDS} $sOptions |" or die "Error: array ID filtering failed: $!\n";
   while (<ARRAYFILT>){
      if ($. == $hCdtPar{sync}{arraytree_id_row} ){
         s/[\n\r]//g;
         my @asLine = split /\t/, $_, -1;
         my $nArrayCount = scalar(@asLine) - $hCdtPar{sync}{array_offset} + 1;
         if ($nArrayCount == scalar(@{$hCdtPar{to}{arraytree_ids}})){
            @asLine[$hCdtPar{sync}{array_offset}-1 .. $#asLine] = @{$hCdtPar{to}{arraytree_ids}};
            $fhTmp->print(join("\t", @asLine), "\n");
         }
         else{
            die "Error: could not sync files because '$sSyncFile' is missing arrayIDs relative to '$sToFile'\n"
         }
      }
      else{
         $fhTmp->print($_);
      }
   }
   close ARRAYFILT;
   close $fhTmp;
   $sStagingFile = $sTmp;
}

# At this point the cdt file was successfully mapped so we can return the result files
my ($sOutName,$sOutPath,$sOutSuffix)    = fileparse($sOut,".cdt");
my ($sSyncName,$sSyncPath,$sSyncSuffix) = fileparse($sSyncFile,".cdt");
my ($sToName,$sToPath,$sToSuffix)       = fileparse($sToFile,".cdt");
copy($sStagingFile, "$sOutPath/$sOutName.cdt");

# Get the accompanying gene tree if available
if ( ($sSyncOn eq "rows") or ($sSyncOn eq "rows-and-columns") ){
   copy("$sToPath/$sToName.gtr", "$sOutPath/$sOutName.gtr") if (-e "$sToPath/$sToName.gtr");
}
else{
   copy("$sSyncPath/$sSyncName.gtr", "$sOutPath/$sOutName.gtr") if (-e "$sSyncPath/$sSyncName.gtr");
}

# Get the accompaning array tree if available
if ( ($sSyncOn eq "columns") or ($sSyncOn eq "rows-and-columns") ){
   copy("$sToPath/$sToName.atr", "$sOutPath/$sOutName.atr") if (-e "$sToPath/$sToName.atr");
}
else{
   copy("$sSyncPath/$sSyncName.atr", "$sOutPath/$sOutName.atr") if (-e "$sSyncPath/$sSyncName.atr");
}

#################
## SUBROUTINES ##
#################

# get_cdt_params
#
# Get line and column IDs for header lines
sub get_cdt_params {
   my ($sFile) = @_;
   open IN, $sFile or die "Error: can't open $sFile: $!\n";
   my %hIDs;
   while (<IN>){
      s/[\n\r]//g;
      my @asLine = split /\t/, $_, -1;
      if (/^GID\tYORF\tNAME\t.*GWEIGHT\t/){
         $hIDs{array_id_row}    = $.;
         $hIDs{gene_id_col}     = 2;
         $hIDs{genetree_id_col} = 1;
         for (my $i=0 ; $i<@asLine ; $i++){
            if ($asLine[$i] eq 'GWEIGHT'){
               $hIDs{array_offset} = $i+1;
            }
         }
      }
      if (/^AID\t.*\t.*\t.*\tARRY\d+X/){
         $hIDs{arraytree_id_row} = $.;
         $hIDs{arraytree_ids}    = [@asLine[$hIDs{array_offset}-1 .. $#asLine]];
      }
      if (/^GENE\d+X\t/){
         $hIDs{gene_offset} ||= $.;
         push @{$hIDs{genetree_ids}}, $asLine[0];
      }
   }
   close IN;
   
   # Check whether all required parameters were set
   my @asKeys;
   if ($sSyncOn eq 'rows'){
      @asKeys = qw(array_id_row array_offset gene_id_col genetree_id_col gene_offset genetree_ids);
   }
   elsif ($sSyncOn eq 'columns'){
      @asKeys = qw(array_id_row arraytree_id_row array_offset arraytree_ids gene_id_col gene_offset);
   }
   elsif ($sSyncOn eq 'rows-and-columns'){
      @asKeys = qw(array_id_row arraytree_id_row array_offset arraytree_ids gene_id_col genetree_id_col gene_offset genetree_ids);
   }
   else{
      die "Error: \n";
   }
   
   foreach my $sKey (@asKeys){
      unless (exists $hIDs{$sKey}){
         die "Error: '$sFile' does not have the required format: missing key '$sKey'\n";
      }
   }
   
   return \%hIDs;
}


# INTERRUPT
#
# Interrupt routine, make sure we exit gracefully for tmp file cleanup
sub INTERRUPT{
   exit(1); # This will call END
}

