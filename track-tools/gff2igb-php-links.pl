#!/usr/bin/env perl

# 04.07.2009 15:42:27 EDT
# Harm van Bakel <hvbakel@gmail.com>
# Converts a gff file into a html file for IGB

# MODULES
use strict;
use warnings;
use Getopt::Long;
use CGI qw(:standard);
use ValidateFiletypes qw(check_gff);

# GET PARAMETERS
my $sHelp      = 0;
my $sInput     = '';
my $sIgbGenome = '';
my $nFlankSize = 100;
GetOptions("help!"    => \$sHelp,
           "input=s"  => \$sInput,
           "genome=s" => \$sIgbGenome,
           "flank:i"  => \$nFlankSize);

# PRINT HELP
$sHelp = 1 unless($sInput and $sIgbGenome);
if ($sHelp) {
   my $sScriptName = ($0 =~ /^.*\/(.+$)/) ? $1 : $0;
   die <<HELP

   Usage: $sScriptName -i <input file> -g <igb-genome-name>
    
    -i <string>
      Input file
    -g <string>
      IGB genome name
    -f <integer>
      Indicates the number of nucleotides to add to either end of each 
      feature to set the display window in IGB.
      default: $nFlankSize
    -help
      This help message
   
HELP
}


###########
## START ##
###########

# Check for proper gff formatting
my @asGffErrors = check_gff($sInput);
if (@asGffErrors){
   unshift @asGffErrors, "Errors in gff file:";
   die join("\n", @asGffErrors), "\n";
}

# Print html header and title
print start_html(-title=>"IGB lookup table for genome '$sIgbGenome'");

# Print the PHP part that checks whether the form is submitted
print "\n\n", '<?php 
if (isset($_POST["SUBMIT"])) {
   echo "<pre>";
   for ($i=0; $i<count($_POST["SELECTION"]);$i++) {
      $id = $_POST["SELECTION"][$i];
      echo $_POST["CHR_$id"] . "\t";
      echo $_POST["SRC_$id"] . "\t";
      echo $_POST["FTR_$id"] . "\t";
      echo $_POST["STA_$id"] . "\t";
      echo $_POST["END_$id"] . "\t";
      echo $_POST["SCR_$id"] . "\t";
      echo $_POST["STR_$id"] . "\t";
      echo $_POST["FRM_$id"] . "\t";
      echo "$id\n";
   }
   echo "</pre>";
}
else{
?>', "\n\n";

# Print the input form
my %hUniqueCheck;
print "<form name='process' method='POST' action='", '<?php echo $_SERVER["PHP_SELF"];?>', "'>\n";
print "<table>";
print "<tr><th align='left'>Selection</th><th align='left'>IGB link</th><th align='left'>Start</th><th align='left'>End</th><th></th></tr>";
open INPUT, $sInput or die "Can't open '$sInput': $!\n";
while (<INPUT>){
   next if /^\s*$/;
   next if /^\s*#/;
   next if /^[Tt]rack/;
   s/[\n\r]$//g;
   my ($sChr, $sSource, $sFeature, $nStart, $nEnd, $nScore, $sStrand, $sFrame, $sID) = split /\t/;

   # Make sure that IDs are unique
   $sID =~ s/ //g;
   die "Error: duplicate rows found for feature'$sID'\n" if (exists $hUniqueCheck{$sID});
   $hUniqueCheck{$sID}++;
   
   # Print table rows
   my $nStartFlank = ($nStart-$nFlankSize)>0 ? $nStart-$nFlankSize : 1;
   my $nEndFlank   = $nEnd+$nFlankSize;
   my $sHref = "http://localhost:7085/UnibrowControl?seqid=$sChr&version=$sIgbGenome&start=$nStartFlank&end=$nEndFlank";
   print "<tr>\n";
   print "    <td><input type='checkbox' name='SELECTION[]' value='$sID'></td>\n";
   print "    <td><a href='$sHref' target='_blank'>$sID</a></td>\n";
   #print "    <td><input type='text' name='STA_$sID' value='$nStart' size='5'></td>\n";
   #print "    <td><input type='text' name='END_$sID' value='$nEnd' size='5'></td>\n";
   print "    <td><input type='hidden' name='CHR_$sID' value='$sChr'>
                  <input type='hidden' name='SRC_$sID' value='$sSource'>
                  <input type='hidden' name='FTR_$sID' value='$sFeature'>
                  <input type='hidden' name='SCR_$sID' value='$nScore'>
                  <input type='hidden' name='STR_$sID' value='$sStrand'>
                  <input type='hidden' name='FRM_$sID' value='$sFrame'>
               </td>\n";
   print "</tr>\n";
}
close INPUT;
print "</table>\n";
print "<input type='submit' value='Update selection' name='SUBMIT'>\n";
print "</form>";

# Closing bracket for the input form part
print "\n<?php } ?>";

print end_html(), "\n";
