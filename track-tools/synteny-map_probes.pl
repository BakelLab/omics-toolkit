#!/usr/bin/env perl


# Script to map tiling array probe measurements from a syntenic genome
# onto another source genome.


# MODULES
use strict;
use Getopt::Long;
use FileHandle;


# GET ARGUMENTS
my $sHelp         = 0;
my $sInputAxt     = '';
GetOptions("help!"        => \$sHelp,
           "input:s"      => \$sInputAxt);


# PRINT HELP
$sHelp = 1 unless($sInputAxt);
if ($sHelp) {
    die <<HELP

    synteny-map_probes.pl -i <axt input file> <probe-file1.txt> ... <probe-fileN.txt>

    Script to map tiling array probe measurements from a syntenic target genome
    onto a reference chromosome.

    arguments:
    -i Input file
        Input file in axt-format containing synteny mappings of the target genome
        against one source chromosome of the reference genome. The axt format is
        described in http://genome.ucsc.edu/goldenPath/help/axt.html
    -h Help

HELP
}

###########
## START ##
###########

# Read synteny mappings
print STDERR "Reading synteny mappings\n";
my ($rhSyntenyMappings, $sRefChr) = &read_axt_mappings($sInputAxt);

# Read all probe measurements
my ($rRemappedProbeMeasurements) = &map_syntenic_measurements(\@ARGV, $rhSyntenyMappings);

# Now sort the mapped positions and write the output
print STDERR "Writing output\n";
&write_output($rRemappedProbeMeasurements, $sRefChr);


#################
## SUBROUTINES ##
#################


#----------------------------------------------------------------------------
# read_axt_mappings($sInputAxt)
#
# Reads the axt-formatted input file with synteny mappings for the 
# reference genome and returns a lookup table for all mouse syntenic regions.

sub read_axt_mappings{
   my ($sInputAxt) = @_;
   my %hSyntenyMappings;
   my $sReferenceChromosome = '';
  
   # Open file
   my $fhAxt = new FileHandle;
   $fhAxt->open($sInputAxt) or die "Can't open axt input file: $!\n";
   while(my $sLine = $fhAxt->getline()){
      next if ($sLine =~ /^\s*$/);
      next if ($sLine =~ /^#/);
      chomp $sLine;
      if ($sLine =~ /^\d+/){
         my ($nBlockID, $sRefChr, $nRefStart, $nRefEnd, $sSyntChr, $nSyntStart, $nSyntEnd, $sStrand, $nScore) = split /\s/, $sLine;
         
         # Check the reference chromosome
         if ($sReferenceChromosome){
            die "Multiple reference chromosomes found in the input file, please process only one chromosome at a time\n" unless ($sReferenceChromosome eq $sRefChr);
         }
         else{
            $sReferenceChromosome = $sRefChr;
         }
        
         # Get the finemapping details
         my $sRefAlignment  = $fhAxt->getline(); chomp $sRefAlignment;
         my $sSyntAlignment = $fhAxt->getline(); chomp $sSyntAlignment;
         my @anFineMappings = &get_finemapping_positions($nRefStart, $sRefAlignment, $sSyntAlignment);
         
         # Store the syntenic region feature for fast hash lookup
         for (my $i=$nSyntStart ; $i<=$nSyntEnd ; $i++){ $hSyntenyMappings{$sSyntChr}{$i} = shift @anFineMappings }
      }
   }
   
   return(\%hSyntenyMappings, $sReferenceChromosome);
}


#----------------------------------------------------------------------------
# get_finemapping_positions($nReferenceStart, $sReferenceAlignment, $sSyntenyAlignment)
#
# Takes the two aligned sequences in the alignment block (the reference and the
# syntenic sequence) and processes each position in the syntenic sequence to
# get the exact mapping to the reference sequence. It returns an array where each
# element ID corresponds to an aligned residue from the syntenic sequence (excluding gaps)
# The value of each array element indicates the corresponding position in the
# reference genome.

sub get_finemapping_positions{
   my ($nRefCounter, $sRefAlignment, $sSyntAlignment) = @_;
   my @anFineMapping;
   
   # Check if the passed strings are of the same length
   die "Found unequal alignment string lengths, please check file format: $!\n" unless (length($sRefAlignment) == length($sSyntAlignment));
   
   for (my $i=0 ; $i<length($sRefAlignment) ; $i++ ){
      push @anFineMapping, $nRefCounter unless (substr($sSyntAlignment, $i, 1) eq '-');
      $nRefCounter++ unless (substr($sRefAlignment, $i, 1) eq '-');
   }
   
   return (@anFineMapping);
}


#----------------------------------------------------------------------------
# map_syntenic_measurements(\@ARGV, $rhSyntenyMappings);
#
# Processes all files with probe measurements (generated with bar2txt) and
# maps probes that fall into syntenic regions with the reference chromosome to their
# corresponding position in that reference chromosome

sub map_syntenic_measurements{
   my ($raProbeInputFiles, $rhSyntenyMappings) = @_;
   my @aaRemappedMeasurements;
   
   my $nCounter = 0;
   print STDERR "Processing probe measurement files, each dot represents 1 million measurements\n";
   foreach my $sFile (@$raProbeInputFiles){
      open INPUT, $sFile or die "Can't open microarray measurement file: $!\n";
      while (<INPUT>){
         next if /^#/;
         next if /^\s*$/;
         chomp;
         my ($sChr, $nPosition, $nMeasurement) = split /\t/;
         if (exists($rhSyntenyMappings->{$sChr}{$nPosition})){
            push @aaRemappedMeasurements, [$rhSyntenyMappings->{$sChr}{$nPosition}, $nMeasurement];
         }
         
         # Print dots so that we can track progress
         if($nCounter==1000000){
            print STDERR '.';
            $nCounter = 0;
         }
         $nCounter++;
         
      }
      close INPUT;
   }
   print STDERR "\n";
   
   return (\@aaRemappedMeasurements);
}



#----------------------------------------------------------------------------
# write_output($rRemappedMeasurements, $sRefChr)
#
# Sort the output by start position and write a gff file

sub write_output{
   my ($rRemappedMeasurements, $sRefChr) = @_;
   
   my @aaSortedMappings = sort { $a->[0] <=> $b->[0] } @$rRemappedMeasurements;
   foreach my $rMeasurement (@aaSortedMappings){
      print join("\t", $sRefChr, @$rMeasurement), "\n";
   }
   
}





