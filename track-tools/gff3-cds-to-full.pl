#!/usr/bin/env perl

# 26.01.2015 16:42:02 EST
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

   Usage: $sScriptName <gff3-file>

   Convert a basic gff3 file containing only CDS annotations (i.e.
   as produced by RAST) to full specification with gene and exon
   descriptions and parent links.
 
   Arguments:
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

my %hIDcheck;
open IN, $ARGV[0] or die "Error: can't open '$ARGV[0]': $!\n";
print "##gff-version\t3\n";
while (<IN>){
   next if (/^\s*$/);
   next if (/^ *#/);
   s/[\n\r]+$//;
   my @asLine = split /\t/;
   my $sType = $asLine[2];
   if ($sType =~ /(CDS|tRNA|rRNA)/){
      my %hAnn = %{gtf_annots_to_hash($asLine[$#asLine])};
      if (exists $hAnn{ID}){
         if (exists $hIDcheck{$hAnn{ID}}){
            warn "Warning: skipping multi-features for '$hAnn{ID}'\n";
         }
         else{
            $hIDcheck{$hAnn{ID}}++;
         
            # Print gene record
            my %hOutAnn = %hAnn;
            $hOutAnn{ID}      = "gene|$hAnn{ID}";
            $asLine[2]        = "gene";
            $asLine[$#asLine] = hash_to_gtf_annots(\%hOutAnn);
            print join("\t", @asLine), "\n";
            
            # Print mRNA record
            %hOutAnn = %hAnn;
            $hOutAnn{ID}      = "mRNA|$hAnn{ID}";
            $hOutAnn{Parent}  = "gene|$hAnn{ID}";
            $asLine[2]        = "mRNA";
            $asLine[$#asLine] = hash_to_gtf_annots(\%hOutAnn);
            print join("\t", @asLine), "\n";
            
            # Print exon record
            %hOutAnn = %hAnn;
            $hOutAnn{ID}      = "exon|$hAnn{ID}";
            $hOutAnn{Parent}  = "mRNA|$hAnn{ID}";
            $asLine[2]        = "exon";
            $asLine[$#asLine] = hash_to_gtf_annots(\%hOutAnn);
            print join("\t", @asLine), "\n";
            
            if ($sType eq "CDS"){
               # Print CDS record
               %hOutAnn = %hAnn;
               $hOutAnn{ID}      = "CDS|$hAnn{ID}";
               $hOutAnn{Parent}  = "mRNA|$hAnn{ID}";
               $asLine[2]        = "CDS";
               $asLine[$#asLine] = hash_to_gtf_annots(\%hOutAnn);
               print join("\t", @asLine), "\n";
            }
         }
      }
      else{
         die "Error: line $. doesn't have an ID annotation\n";
      }
   }
}
close IN;



#################
## SUBROUTINES ##
#################

# gtf_annots_to_hash
#
# Parse gtf annotations and return key-value pairs
sub gtf_annots_to_hash {
   my ($sAnnots) = @_;
   my %hReturn;
   
   my @asPairs = split /;/, $sAnnots;
   foreach my $sPair (@asPairs){
      my ($sKey, $sVal) = split /=/, $sPair;
      if ($sKey eq 'tag'){
         push @{$hReturn{$sKey}}, $sVal;
      }
      else{
         die "Error: Duplicate key entry '$sKey' found with values '$hReturn{$sKey}' and '$sVal'\n" if (exists $hReturn{$sKey});
         $hReturn{$sKey} = $sVal;
      }
   }
   return \%hReturn;
}

# hash_to_gtf_annots
#
# hash_to_gtf_annots
sub hash_to_gtf_annots {
   my ($rAnn) = @_;
   my @asOut;
   foreach my $sKey (sort keys %$rAnn){
      push @asOut, $sKey . "=" . $rAnn->{$sKey};
   }
   return join(";", @asOut);
}
