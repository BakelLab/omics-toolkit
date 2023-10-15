#!/usr/bin/env perl

# 09.03.2011 14:01:57 EST
# Harm van Bakel <hvbakel@gmail.com>

# MODULES
use strict;
use warnings;
use Getopt::Long;
use File::Temp qw(tempfile tempdir);

# GLOBALS
$ENV{GETORF} ||= 'getorf';
$ENV{TMPDIR} ||= "/tmp";
$SIG{'INT'}=$SIG{'HUP'}=$SIG{'ABRT'}=$SIG{'QUIT'}=$SIG{'TRAP'}=$SIG{'STOP'}=\&INTERRUPT;

# GET PARAMETERS
my $sHelp        = 0;
my $nMinLength   = 45;
my $flIncludeORF = 0;
GetOptions("help!"     => \$sHelp,
           "minorf:i"  => \$nMinLength,
           "orf!"      => \$flIncludeORF);
my $sInput = shift @ARGV;

# PRINT HELP
$sHelp = 1 unless($sInput);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName [options] <fasta-file>
   
   When provided with a multifasta file, this script
   will return the length of each sequence, and optionally
   the length of the largest ORF on the forward and reverse
   strand.
    
   Options:
    -o --orf
      Report largest ORF on forward and reverse strand.
    -m --minorf <integer>
      Minimum length of reported ORFs. Default: $nMinLength
    -help
      This help message
   
HELP
}


##########
## MAIN ##
##########

# Reformat the fasta input for easier getORF parsing
my ($sTmpIn, $raSeqLengths) = convert_fasta($sInput);


unless ($flIncludeORF){
   print "#id\tlength\tN-count\tGC-content\tN-content\n";
   foreach my $rSeq (@$raSeqLengths){
      my ($sID, $nLength, $nNcount, $nGCcontent, $nNcontent) = @$rSeq;
      $sID =~ s/^>//;
      $sID =~ s/___/;/g;
      my $sGCcontent = $nGCcontent eq 'NA' ? 'NA' : sprintf "%.5f", $nGCcontent;
      print "$sID\t$nLength\t$nNcount\t$nGCcontent\t$nNcontent\n";
   }
}
else{
   # Run getORF and parse output into an array
   my $raGetORF = run_getORF($sTmpIn, $nMinLength);
   
   my %hORFs;
   foreach my $rRecord (@$raGetORF){
      my ($sHeader, $sSeq) = @$rRecord;
      my ($sID, $nOstart, $nOend, $nOrientation, $nTlength) = ('',0,0,0,0);
      
      # Check header format
      if ($sHeader =~ /^(.*) \[(\d+) \- (\d+)\] (\d+)$/){
         ($sID, $nOstart, $nOend, $nOrientation, $nTlength) = ($1,$2,$3,0,$4);
      }
      elsif ($sHeader =~ /^(.*) \[(\d+) \- (\d+)\] \(REVERSE SENSE\) (\d+)$/){
         ($sID, $nOstart, $nOend, $nOrientation, $nTlength) = ($1,$3,$2,1,$4);
      }
      else{
         die "Error: incompatible getORF header: '$_'\n";
      }
      
      # Reformat some stuff
      $sID =~ s/^>//;
      $sID =~ s/_\d+$//;
      $sID =~ s/___/;/g;
      my $nOlength = $nOend - $nOstart + 1;
      
      if (exists $hORFs{$sID}{$nOrientation}){
         if ($nOlength > $hORFs{$sID}{$nOrientation}{length}){
            $hORFs{$sID}{$nOrientation}{length} = $nOlength;
            $hORFs{$sID}{$nOrientation}{seq}    = $sSeq;
         }
         
      }
      else{
         $hORFs{$sID}{$nOrientation}{length} = $nOlength;
         $hORFs{$sID}{$nOrientation}{seq}    = $sSeq;
      }
   }
   close GETORF;
   
   
   # Now have a final run through the data, making sure to only print
   # orientation info if there is a clear ORF bias towards one of the strands
   foreach my $rRecord (@$raSeqLengths){
      my ($sID, $nSeqLength, $nNcount, $nGCcontent, $nNcontent) = @$rRecord;
      my $nORFfwdLength = exists($hORFs{$sID}{0}) ? $hORFs{$sID}{0}{length} : 0;
      my $nORFrevLength = exists($hORFs{$sID}{1}) ? $hORFs{$sID}{1}{length} : 0;
      my $sORFfwdSeq    = exists($hORFs{$sID}{0}) ? $hORFs{$sID}{0}{seq} : "";
      my $sORFrevSeq    = exists($hORFs{$sID}{1}) ? $hORFs{$sID}{1}{seq} : "";
      
      # Final check to make sure everything is ok
      die "Error: FWD orf length and sequence length don't match: $nORFfwdLength\t$sORFfwdSeq\n" unless (3*(length($sORFfwdSeq)) == $nORFfwdLength);
      die "Error: REV orf length and sequence length don't match: $nORFrevLength\t$sORFrevSeq\n" unless (3*(length($sORFrevSeq)) == $nORFrevLength);
      
      # Print output
      print join("\t", $sID, $nSeqLength, $nNcount, $nGCcontent, $nNcontent, $nORFfwdLength, $nORFrevLength, $sORFfwdSeq,$sORFrevSeq), "\n";
   }
}

#################
## SUBROUTINES ##
#################

# convert_fasta
#
# Rewrite the fasta file for getorf
sub convert_fasta {
   my ($sFastaIn) = @_;
   my @aaSeqLengths;

   # Process the INPUT
   my $sFastaHeader = '';
   my $sFastaSeq    = '';
   my ($fhTmp, $sTmp) = tempfile('getORF-XXXXX', DIR=>$ENV{TMPDIR}, UNLINK=>1);
   open INPUT, "<$sFastaIn" or die "Error: can't read the fasta file\n";
   while (<INPUT>){
      s/[\n\r]+$//;
      if (/^>/ or eof){
         if (eof){
            die "Error: file ends in fasta header without sequence\n" if (/^>/);
            $sFastaSeq .= $_;
         }
         $sFastaSeq  =~ s/\s//g;
         if ($sFastaHeader){
            my $nSeqLength = length($sFastaSeq);
            $sFastaSeq =~ s/.{100}/$&\n/sg;
            $sFastaSeq =~ s/\n$//; # Make sure to remove any trailing newlines
            print $fhTmp join("\n", "$sFastaHeader $nSeqLength", $sFastaSeq), "\n";

            # Keep fasta lengths
            my $sID = $sFastaHeader;
            $sID =~ s/^>//;
            $sID =~ s/___/;/g;
            my $nGCcount   = $sFastaSeq =~ tr/GCgc/GCgc/;
            my $nATcount   = $sFastaSeq =~ tr/ATat/ATat/;
            my $nNcount    = $sFastaSeq =~ tr/Nn/Nn/;
            my $nGCcontent = ($nGCcount+$nATcount)>0 ? $nGCcount/($nGCcount+$nATcount) : "NA";
            my $nNcontent  = ($nNcount)>0 ? $nNcount/($nGCcount+$nATcount+$nNcount) : "NA";
            push @aaSeqLengths, ([$sID, $nSeqLength, $nNcount, $nGCcontent, $nNcontent]);
         }
         $sFastaHeader = $_;
         $sFastaHeader =~ s/\s+.*$//;
         $sFastaHeader =~ s/;/___/g;
         $sFastaSeq    = "";
      }
      else{
         next if (/^\s*$/);
         next if (/^ *#/);
         $sFastaSeq .= $_ if ($sFastaHeader);
      }
   }
   close INPUT;
   $fhTmp->close();
   return $sTmp, \@aaSeqLengths;
}


# run_getORF
#
# Run getORF and get the output in an array
sub run_getORF {
   my ($sInput, $nMinLength) = @_;
   
   my @asReturn;
   my $sFastaHeader = '';
   my $sFastaSeq    = '';
   open GETORF, "$ENV{GETORF} -minsize $nMinLength -sequence $sInput -stdout -auto |" or die "Error: can't run getorf: $!\n";
   while (<GETORF>){
      s/[\n\r]+$//;
      if (/^>/ or eof){
         $sFastaSeq .= $_ if (eof);
         if ($sFastaHeader){
            push @asReturn, ([$sFastaHeader, $sFastaSeq]);
         }
         $sFastaHeader = $_;
         $sFastaSeq    = "";
      }
      else{
         next if (/^\s*$/);
         $sFastaSeq .= $_ if ($sFastaHeader);
      }
   }
   close INPUT;
   return \@asReturn;
}


# INTERRUPT
#
# Interrupt routine, make sure we exit gracefully for tmp file cleanup
sub INTERRUPT{
   exit(1); # This will call END
}
