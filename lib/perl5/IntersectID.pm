#!/usr/bin/perl

# Package with some common functions to read IDs from columns and rows

package IntersectID;

use strict;
use warnings;
use Carp qw( croak );
require Exporter;

our @ISA       = qw(Exporter);
our @EXPORT_OK = qw(read_ids_from_row read_ids_from_col remap_hash_keys_from_file);


# read_ids_from_row
#
# Read unique identifiers from an input file row
# Returns a hash reference with IDs as key and order in which the IDs were 
# first encountered and a hash reference with duplicate IDs and their occurrences
sub read_ids_from_row{
   my %args = (file          => "",
               id_row        => 1,
               offset        => 1,
               fieldsep      => "\t",
               split_char    => "",
               preserve_case => 1,
               @_);

   $args{offset}--; # zero-base the offset to be used as an array index
   my %hLUT;
   my %hDuplicateIDs;
   my $nCount = 0;
   my $nOrderedIDs = 0;
   open LUT, $args{file} or croak "Can't open $args{file}: $!";
   while(<LUT>){
      if ($args{id_row} == $.){
         s/[\n\r]//g;
         my (@asLine) = split /$args{fieldsep}/, $_, -1;
         if ($args{offset} <= $#asLine){
            for( my $i=$args{offset} ; $i<@asLine; $i++){
               my @asIDs;
               if ($args{split_char}){
                  push @asIDs, (split /$args{split_char}/, $asLine[$i]);
               }
               else{
                  push @asIDs, $asLine[$i];
               }
               _append_ids_to_lut(\%hLUT, \%hDuplicateIDs, \$nCount, \@asIDs, $args{preserve_case});
            }
            last;
         }
         else{
            croak "Error: ID file offset is larger than the number of ID fields";
         }
      }
   }
   close LUT;   
   return (\%hLUT, \%hDuplicateIDs);
}


# read_ids_from_col
#
# Read unique identifiers from an input file column
# Returns a hash reference with IDs as key and order in which the IDs were 
# first encountered and a hash reference with duplicate IDs and their occurrences
sub read_ids_from_col{
   my %args = (file          => "",
               id_col        => 1,
               offset        => 1,
               fieldsep      => "\t",
               split_char    => "",
               preserve_case => 1,
               @_);
   
   $args{id_col}--; # zero-base the id col ID to be used as an array index
   my %hLUT;
   my %hDuplicateIDs;
   my $nCount = 0;
   open LUT, $args{file} or croak "Can't open $args{file}: $!";
   while(<LUT>){
      next if ($. < $args{offset});
      next if /^\s*$/;
      next if /^\s*#/;
      s/[\n\r]//g;
      my (@asLine) = split /$args{fieldsep}/, $_, -1;
      if ($args{id_col} < @asLine){
         my $sIDfield = $asLine[$args{id_col}];
         my @asIDs;
         if ($args{split_char}){
            push @asIDs, (split /$args{split_char}/, $sIDfield);
         }
         else{
            push @asIDs, $sIDfield;
         }
         _append_ids_to_lut(\%hLUT, \%hDuplicateIDs, \$nCount, \@asIDs, $args{preserve_case});
      }
      else{
         croak "Error: number of fields in the ID file is lower than the supplied ID column number on line $.";
      }
   }
   close LUT;
   return (\%hLUT, \%hDuplicateIDs);
}


# _append_ids_to_lut
#
# Append IDs to LUT
sub _append_ids_to_lut{
   my ($rLUT, $rDuplicates, $rCount, $rArray, $nPreserveCase) = @_;
   foreach my $sID (@$rArray){
      $sID =~ s/^\s+//; # Remove leading whitespace
      $sID =~ s/\s+$//; # Remove trailing whitespace
      $sID = lc($sID) unless ($nPreserveCase);
      if (exists $rLUT->{$sID}){
         $rDuplicates->{$sID}++;
      }
      else{
         $rLUT->{$sID} = $$rCount;
      }
      $$rCount++;
   }
}


# remap_hash_keys_from_file
#
# Remap the lut according to a file with id mappings between the
# ID and filter file
sub remap_hash_keys_from_file {
   my %args = (lut           => undef,
               file          => "",
               offset        => 1,
               preserve_case => 1,
               fieldsep      => "\t",
               reverse_map   => 0,
               @_);
   
   # Read the mappings into a LUT
   my %hMappingLUT;
   my %hUniqueMapValCheck;
   open MAPPING, $args{file} or croak "Error: can't open mapping file '$args{file}': $!";
   while (<MAPPING>){
      next if ($. < $args{offset});
      next if (/^\s*$/);
      next if (/^\s*#/);
      s/[\n\r]+$//;
      my $sLine = $args{preserve_case} ? $_ : lc($_);
      my (@asLine) = split /$args{fieldsep}/, $sLine, -1;
      croak "Error: mapping file does not have a paired ID mapping entry on line $." unless (@asLine == 2);
      my ($sMapKey, $sMapVal) = @asLine;
      ($sMapKey, $sMapVal)    = ($sMapVal, $sMapKey) if ($args{reverse_map});
      if (exists $hMappingLUT{$sMapKey}){
         croak "Error: duplicate keys found in mapping file, line $.";
      }
      elsif (exists $hUniqueMapValCheck{$sMapVal}){
         croak "Error: duplicate values found in mapping file, line $.";
      }
      else{
         $hMappingLUT{$sMapKey} = $sMapVal;
      }
   }
   close MAPPING;
   
   # Remap the LUT keys to the new key values
   my $nRemapCount = 0;
   my %hNewLUT;
   foreach my $sKey ( keys(%{$args{lut}}) ){
      my $sCompKey = $args{preserve_case} ? $sKey : lc($sKey);
      if (exists $hMappingLUT{$sCompKey}){
         $hNewLUT{$hMappingLUT{$sCompKey}} = $args{lut}{$sKey};
         $nRemapCount++;
      }
      else{
         $hNewLUT{$sKey} = $args{lut}{$sKey};
      }
   }
   my $nTotalCount = scalar(keys(%hNewLUT));
   if ($nTotalCount == $nRemapCount){
      warn "Remapped all IDs according to mapping file\n";
   }
   else{
      warn "Remapped $nRemapCount out of $nTotalCount IDs according to mapping file\n";
   }
   return(\%hNewLUT);
}


1;
