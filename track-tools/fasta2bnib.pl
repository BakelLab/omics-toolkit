#!/usr/bin/env perl

# Convert fasta files to binary sequence files (bnib) that can be loaded into the IGB

# MODULES
use strict;
use FileHandle;

# GLOBALS
$ENV{JAVA_BIN}      ||= '/usr/bin/java';
$ENV{IGB_JAR_PATH}  ||= '/home/projects/hugheslab/opt/igb-jarfiles';
$ENV{BNIB_COMMAND}  ||= join('',
			     $ENV{JAVA_BIN}, ' -classpath ',
			     $ENV{IGB_JAR_PATH}, '/RMgenometry.jar:',
			     $ENV{IGB_JAR_PATH}, '/RMgenoviz.jar:',
			     $ENV{IGB_JAR_PATH}, '/RMigb.jar ',
			     'com.affymetrix.igb.parsers.NibbleResiduesParser ');

# PRINT USAGE
unless (@ARGV) {
    die <<HELP

    Usage: fasta2bnib <fasta-file1> ... <fasta-fileN>

    Converts fasta-formatted chromosome sequence files into binary sequence files for
    loading into IGB. Make sure that the fasta header matches the chromosome names in
    the chromosome info and gff files.
    
HELP
}


foreach my $sFastaFile (@ARGV){
   
   # Read fasta header
   my $fhFasta = new FileHandle;
   $fhFasta->open($sFastaFile) or die "Can't open $sFastaFile: $!\n";
   my $sHeader = $fhFasta->getline();
   $fhFasta->close();
   
   # Check header
   die "Incorrect file format, need fasta-formatted header line\n" unless ($sHeader =~ /^>/);
   $sHeader =~ s/[>\n\r\s\t]//g;
   
   
   my $sBnibName = join('.', $sHeader, 'bnib');
   `$ENV{BNIB_COMMAND} $sHeader $sFastaFile $sBnibName 1`;
}

