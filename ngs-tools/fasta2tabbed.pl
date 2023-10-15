#!/usr/bin/env perl

# 08.08.2011 13:44:05 EDT
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;

# GET PARAMETERS
my $sHelp        = 0;
GetOptions("help!"   => \$sHelp);

# PRINT HELP
$sHelp = 1 unless(@ARGV>0);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName <fasta-file>
   
   Converts a fasta-formatted file to a tab delimited format
   with the header and sequence in separate columns.
    
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

my $raFasta = read_fasta($ARGV[0]);
foreach my $rFasta (@$raFasta){
   print join("\t", @$rFasta), "\n";
}

#################
## SUBROUTINES ##
#################

# read_fasta
#
# Reads fasta file into a hash
sub read_fasta {
   my ($sFASTA) = @_;
   
   my @aaFasta;
   my $sFastaHeader = '';
   my $sFastaSeq    = '';
   open FASTA, "<$sFASTA" or die "Error: can't read the fasta file\n";
   while (<FASTA>){
      s/[\n\r]+$//;
      if (/^>/){
         die "Error: file ends in fasta header without sequence\n" if (eof);
         $sFastaSeq  =~ s/\s//g;
         push @aaFasta, ([$sFastaHeader,$sFastaSeq]) if ($sFastaHeader);
         
         # Reset for the next sequence
         $sFastaHeader = $_;
         $sFastaHeader =~ s/\s*$//;
         $sFastaHeader =~ s/^>\s*//;
         $sFastaSeq    = "";
      }
      elsif (eof){
         $sFastaSeq .= $_;
         $sFastaSeq  =~ s/\s//g;
         push @aaFasta, ([$sFastaHeader,$sFastaSeq]) if ($sFastaHeader);
      }
      else{
         next if (/^\s*$/);
         next if (/^ *#/);
         $sFastaSeq .= $_ if ($sFastaHeader);
      }
   }
   close FASTA;
   return \@aaFasta;
}

