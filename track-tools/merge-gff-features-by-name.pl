#!/usr/bin/env perl

# Merge GFF ORFs/exons into transcripts, with the intron parts included

# MODULES
use strict;
use warnings;

unless(@ARGV){
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die "\n Usage: $sScriptName <gff-file>\n\n";
}
my $sFile = shift @ARGV;


# Open and process entries into a hash of arrays
my @asHeader = qw(chr source type start end score strand frame id);
my %hFeatures;
open IN, $sFile or die "Error: can't open '$sFile': $!\n";
while (<IN>){
   next if /^track/;
   next if /^\s*$/; 
   next if /^#.*$/;
   s/[\n\r]//g;
   my (@asLine) = split /\t/;
   die "Error: incorrect number of fields on line $.\n" unless (@asLine == 9);
   my %hLine;
   @hLine{@asHeader} = @asLine;
   
   die "Error: start position needs to be an integer on line $.\n"  unless ($hLine{start}  =~ /\d+/);
   die "Error: end position needs to be an integer on line $.\n"    unless ($hLine{end}    =~ /\d+/);
   die "Error: unknown strand type on line $.\n"                    unless ($hLine{strand} =~ /[+-.]/);
   
   if (exists($hFeatures{$hLine{id}})){
      die "Error: Feature $hLine{id} was found on two strands\n"       unless ($hFeatures{$hLine{id}}{strand} eq $hLine{strand});
      die "Error: Feature $hLine{id} was found on more than one chr\n" unless ($hFeatures{$hLine{id}}{chr} eq $hLine{chr});
      $hFeatures{$hLine{id}}{start} = $hLine{start} if ($hLine{start} < $hFeatures{$hLine{id}}{start});
      $hFeatures{$hLine{id}}{end}   = $hLine{end}   if ($hLine{end}   > $hFeatures{$hLine{id}}{end});
      $hFeatures{$hLine{id}}{score} = $hLine{score} if ($hLine{score} > $hFeatures{$hLine{id}}{score});
   }
   else{
      $hFeatures{$hLine{id}} = {%hLine};
   }
}
close IN;

# Print the merged file, after sorting by chromosome and start site
my @asSorted = sort {$hFeatures{$a}{chr} cmp $hFeatures{$b}{chr} || $hFeatures{$a}{start} <=> $hFeatures{$b}{start} } keys(%hFeatures);
foreach my $sKey (@asSorted){
   print join("\t", @{$hFeatures{$sKey}}{@asHeader}), "\n";
}
