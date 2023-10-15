#!/usr/bin/perl

# Basic package to validate various standard biological filetypes
# See individual subroutines to check which checks are supported

package ValidateFiletypes;

use strict;
use warnings;
require Exporter;
our @ISA       = qw(Exporter);
our @EXPORT_OK = qw(check_bed check_gff);


# check_bed
#
# Make sure the input file conforms to the bed format specifications
sub check_bed {
   my $sFile      = shift @_;
   my $nMinFields = 4;
   $nMinFields    = shift @_ if(@_);
   my $nMaxErrors = 15;
   $nMaxErrors    = shift @_ if(@_);
   my @asErrors;
   
   # Check bed file format
   # Warn if an uncommented track line is found since fjoin won't run
   my $nFieldCount = 0;
   open IN, $sFile or die "Can't open '$sFile': $!\n";
   while (<IN>){
      next if /^\s*$/;
      next if /^#.*$/;
      next if (/^[Tt]rack/);  # Note: track name checking isn't supported yet
      s/[\n\r]//g;
      my @asLine = split /\t/;
      
      # Check for the correct number of fields
      push @asErrors, "  Insufficient number of fields on line $.: missing mandatory fields" unless (@asLine >=$nMinFields);
      push @asErrors, "  Insufficient number of fields on line $.: thick end missing"        if (@asLine == 7);
      if ($nFieldCount){
         push @asErrors, "  Field count changes on line $." unless ($nFieldCount == @asLine);
      }
      $nFieldCount = scalar(@asLine) if (scalar(@asLine)>$nFieldCount);
      
      # Start checking mandatory fields
      push @asErrors, "  Start site must be a number on line $."    unless($asLine[1] =~ /^\d+$/);
      push @asErrors, "  End site must be a number on line $."      unless($asLine[2] =~ /^\d+$/);
      push @asErrors, "  Start must be smaller than end on line $." unless($asLine[1] <= $asLine[2]);
      
      # Check optional fields
      if (@asLine >= 6){
         push @asErrors, "  Unknown strand type '$asLine[5]' on line $." unless($asLine[5] =~ /^[+-.]$/);
      }
      if (@asLine >= 8){
         push @asErrors, "  Thick start must be a number on line $."               unless($asLine[6] =~ /^\d+$/);
         push @asErrors, "  Thick end site must be a number on line $."            unless($asLine[7] =~ /^\d+$/);
         push @asErrors, "  Thick start must be smaller than thick end on line $." unless($asLine[6] <= $asLine[7]);
      }
      if (@asLine >= 9){
         push @asErrors, "  Incorrect color specification on line $."   unless($asLine[8] =~ /^(\d{1,3},\d{1,3},\d{1,3})|0$/);
      }
      if (@asLine >9){
         if (@asLine >= 12){
            $asLine[10] =~ s/,$//;
            $asLine[11] =~ s/,$//;
            my @anBlockSizes  = split /,/, $asLine[10];
            my @anBlockStarts = split /,/, $asLine[11];
	    my $nDerivedEnd   = $asLine[1] + $anBlockSizes[$#anBlockSizes] + $anBlockStarts[$#anBlockStarts];
            push @asErrors, "  Inconsistency between number of blocks on line $." unless ((@anBlockSizes == $asLine[9]) and (@anBlockStarts == $asLine[9]));
	    push @asErrors, "  'Chr end' does not match 'Chr start + BlockStarts[last] + BlockSizes[last]' on line $." unless($asLine[2] == $nDerivedEnd);
         }
         else{
            push @asErrors, "  Incorrect field count for full bed specification on line $.";
         }
      }
      if(@asErrors > $nMaxErrors){
         $asErrors[$#asErrors] = "--- Further error messages truncated ---";
         last;
      }
   }
   close IN;
   return @asErrors;
}

# check_gff
# 
# Make sure the input file conforms to the gff format specifications
sub check_gff {
   my $sFile      = shift @_;
   my $nMaxErrors = 15;
   $nMaxErrors    = shift @_ if(@_);
   my @asErrors;
   
   # Check gff file format
   # Warn if an uncommented track line is found since fjoin won't run otherwise
   my %hFeatures;
   open IN, $sFile or die "Can't open '$sFile': $!\n";
   while (<IN>){
      next if /^\s*$/;
      next if /^#.*$/;
      next if (/^[Tt]rack/);  # Note: track name checking isn't supported yet
      s/[\n\r]//g;
      my @asLine = split /\t/;
      
      # Check for the correct number of fields
      push @asErrors, "  Insufficient number of fields on line $."        unless (@asLine ==9);
      push @asErrors, "  Start position needs to be a number on line $."  unless ($asLine[3] =~ /^\d+$/);
      push @asErrors, "  End position needs to be a number on line $."    unless ($asLine[4] =~ /^\d+$/);
      push @asErrors, "  Score needs to be a number or '.' on line $."    unless ($asLine[5] =~ /^\d+|\.$/);
      push @asErrors, "  Unknown strand type '$asLine[6]' on line $."     unless ($asLine[6] =~ /^[+-.]$/);
      push @asErrors, "  Frame needs to be between 0-2 or '.' on line $." unless ($asLine[7] =~ /^[012.]$/);
      push @asErrors, "  Start needs to be smaller than end on line $."   unless ($asLine[3] <= $asLine[4]);
      
      # Make sure that features with same name are restricted to same strand and chromosome
      if (exists($hFeatures{$asLine[8]})){
         push @asErrors, "  Feature '$asLine[8]' was found on multiple chromosomes" unless ($asLine[0] eq $hFeatures{$asLine[8]}{chr});
         push @asErrors, "  Feature '$asLine[8]' was found on both strands"         unless ($asLine[6] eq $hFeatures{$asLine[8]}{strand});
      }
      else{
         $hFeatures{$asLine[8]}{chr}    = $asLine[0];
         $hFeatures{$asLine[8]}{strand} = $asLine[6];
      }
      
      if(@asErrors > $nMaxErrors){
         $asErrors[$#asErrors] = "--- Further error messages truncated ---";
         last;
      }
   }
   close IN;
   
   return @asErrors;
}


1;
