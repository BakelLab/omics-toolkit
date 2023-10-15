#!/usr/bin/env perl

# create-igb-genome.pl
# Creates a web-based genome for the Integrated Genome Browser

# MODULES

use strict;
use Getopt::Long;
use FileHandle;
use DBI;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::TranscriptAdaptor;
use Cwd;
use File::Basename;
use File::Copy;
require "igb-conf.pm";

# GET ARGUMENTS
my $sHelp     = 0;
my $sList     = 0;
my @anAnnotIDs;
GetOptions("help!"    => \$sHelp,
	   "list!"    => \$sList,
	   "annot:s"  => \@anAnnotIDs);

# PRINT HELP
$sHelp = 1 unless($sList or scalar(@anAnnotIDs));
if ($sHelp) {
    die <<HELP
    create-igb-genome.pl
    Creates all necessary files for an igb quickload genome
    
    arguments:
    -l list OffBase genomes & reporter libraries
        Gives a list of currently available genomes and reporter libraries in base
    -a Annotation ID
        The reporterlibraryannotation identifier.
    -h Help
HELP
}


######################################
#               START                #
######################################

my $dbhBase = &base_connect();

if ($sList){
    &list_base_genomes($dbhBase);
}

if (scalar(@anAnnotIDs)){
    my ($sOrganism, $sEnsBuild, $sEnsDate, $sOutDir, $sOutDirSgr, $fhLog, $rAnnotMap, $sGenespringMappings) = &check_and_prepare_output($dbhBase, @anAnnotIDs);
    my $reporters = &retrieve_reporters($dbhBase, $fhLog, $sOrganism, @anAnnotIDs);
    my (@asOligoAnnots) = &create_oligo_annotation_files($dbhBase, $fhLog, $sOutDir, $rAnnotMap, @anAnnotIDs);
    my (@asGeneAnnots)  = &create_genome_seq_and_annot_files($fhLog, $sOrganism, $sEnsBuild, $sOutDir);

    # Additional optional track information
    my @asAdditional;
    if ($sOrganism eq 'yeast'){
	my (@asSGDannot)         = &sdg_reannotate_and_peakcluster($fhLog, $sOrganism, $sOutDir, $ENV{GENE_ANNOTS_FILE}, $sGenespringMappings);
	#my (@asSGDtransfac)      = &retrieve_sgd_regulatory_mappings($fhLog, $sOrganism, $sOutDir);
	my (@asHarbisonMappings) = &map_harbison_transfacs($fhLog, $sOrganism, $sOutDir);
	my (@asDavidMappings)    = &map_genome_wide_transcription($fhLog, $sOrganism, $sOutDir);
	my (@asLexAMappings)     = &copy_lexA_control_plasmid_mappings($fhLog, $sOrganism, $sOutDir);
	my (@asRickYoungSignaling) = &map_rickyoung_signaling($fhLog, $sOrganism, $sOutDirSgr, $sEnsDate, $reporters, @anAnnotIDs);
	@asAdditional = (@asSGDannot, @asHarbisonMappings, @asDavidMappings, @asLexAMappings);
    }
    &create_annots_file($fhLog, $sOutDir, @asOligoAnnots, @asGeneAnnots, @asAdditional);
    my (@asDirName) = split /\//, $sOutDir;
    &append_base_quickload_files($sOrganism, $fhLog, $asDirName[$#asDirName]);
    $fhLog->close();
}

# Cleanup
$dbhBase->disconnect();





######################################
#            SUBROUTINES             #
######################################


# base_connect($dbh)
#
# Connect to the base database
# Returns the database handle
sub base_connect {
    my $db  = join('', 'dbi:Pg:host=', $ENV{BASE_HOST}, ';dbname=', $ENV{BASE_DB});
    my $dbh = DBI->connect($db, $ENV{BASE_USER}, $ENV{BASE_PASSWD} );
    die "Failed to connect to database\n" if !$dbh;
    return($dbh);
}


# list_base_genomes($dbh)
#
# Print out a listing of the genomes available in base
sub list_base_genomes {
    my $dbh = shift @_;
    
    # Get listing
    my $sSQL          = join (' ',
				    'SELECT A."id", R."name", A."addedDate", A."name"',
				    'FROM "reporterlibraryannotation" A',
				    'LEFT JOIN "reporterlibrary" R ON A."library"=R."id"',
				    'ORDER BY R."name", A."name"');
    my $oQuery        = $dbh->prepare($sSQL);
    my $n             = $oQuery->execute();
    
    if ( $n eq '0E0' ){
	die "Nothing found \n";
    }
    else{
	printf "%-5s  %-30s %-12s %-10s", 'ID', 'NAME', 'DATE', 'ENSEMBL_BUILD';	
	print "\n------------------------------------------------------------------------\n";
	while (my @asRow = $oQuery->fetchrow_array() ){
	    printf ("%-5s  %-30s %-12s %-10s", @asRow);
	    print "\n";
	}
	print "\n\nUse the IDs from the first colum to generate the appropriate genome.\n";
    }
    
    #Cleanup
    $oQuery->finish();
}



# check_and_prepare_output ($dbh, @annotIDs)
# 
# Validates the data sources and the inputted data
# Returns various needed info
sub check_and_prepare_output{
    my ($dbh, @anAnnotations) = @_;
    my $sOrganism;
    my $sEnsBuild;
    my $sEnsDate;
    my %hAnnotMap;
    
    # Check the input (only annotations from 1 species and 1 Ensembl build; species config)
    my $sAnnotations  = join ('', "'", join ("','", @anAnnotations), "'");
    my $sSQL          = join ('',
			      'SELECT R."name", A."name", A."id", A."assemblyDate" ',
			      'FROM "reporterlibraryannotation" A ',
			      'LEFT JOIN "reporterlibrary" R ON A."library"=R."id" ',
			      'WHERE A."id" IN (', $sAnnotations,') ');
    my $oQuery        = $dbh->prepare($sSQL);
    my $n             = $oQuery->execute();
    if ( $n eq '0E0' ){
	die "Nothing found \n";
    }
    else{
	die "Error: Could not find all supplied annotations in base!\n" if ($n < scalar(@anAnnotations) );
	while (my @asRow = $oQuery->fetchrow_array() ){
	    my ($sLibName, $sTmpBuild, $nAssemblyID, $sTmpDate) = @asRow;
	    my ($tmp, $sTmpOrg, @asRest) = split /\s/, $sLibName;
	    $sOrganism = $sTmpOrg   unless ($sOrganism);
	    $sEnsBuild = $sTmpBuild unless ($sEnsBuild);
	    $sEnsDate  = $sTmpDate unless ($sEnsDate);
	    die "Error: Annotations from different organisms selected!\n"      if ($sOrganism ne $sTmpOrg);
	    die "Error: Annotations from different Ensembl builds selected!\n" if ($sEnsBuild ne $sTmpBuild);
	    die "Error: Annotations from different Ensembl build dates selected!\n" if ($sEnsDate ne $sTmpDate);
	    $hAnnotMap{$nAssemblyID} = $sLibName;
	}
    }
    die "Error: No entry found in config file for organism '$sOrganism'!\n" unless(exists($ENV{ORGANISM}{$sOrganism}));
    $oQuery->finish();

    # Verify if the output dir already exists, don't overwrite previous genomes!
    my $sOutName = join('_', $sOrganism, $sEnsDate, $sEnsBuild);
    my $sOutDir  = join('/', cwd(), $sOutName);
    if(-e "$sOutDir"){ die "Output folder '$sOutDir' already exists, please remove first!"; }
    else             { mkdir("$sOutDir") or die "Could not create output dir $sOutDir: $!\n"};

    # create subdir for .sgr files
    my $sOutDirSgr = $sOutDir . "/sgrfiles";
    mkdir("$sOutDirSgr") or die "Could not create sgr output dir $sOutDirSgr: $!\n";
    
    # Check if a genespring mapping file exists for this genome
    my $sGenespringMappings = "$sOutName.mappings";
    die "No genespring mappings file found for genome $sOutName!" unless (-e "$ENV{GENESPRING_MAPPINGS_DIR}/$sGenespringMappings");
    
    # Open a logfile with all script output
    my $sOligoAnnots = join('; ', values(%hAnnotMap));
    my $sDate        = `date`;
    my $fhLog        = new FileHandle;
    $fhLog->open(">$sOutDir/$ENV{LOG_FILENAME}") or die "Could not open $ENV{LOG_FILENAME}: $!\n";
    $fhLog->print(join("\n", "Creation date:\t$sDate", "Organism:\t$sOrganism", "Ensembl build:\t$sEnsBuild",
		       "Ensembl build date:\t$sEnsDate", "Oligo tracks:\t $sOligoAnnots\n"));

    return($sOrganism, $sEnsBuild, $sEnsDate, $sOutDir, $sOutDirSgr, $fhLog, \%hAnnotMap, $sGenespringMappings);
}


# create_oligo_annotation_files ($dbh, @anAnnotIDs)
#
# Create BED annotation files for the selected reporter library mappings
# Returns a list of the annotation files that it created
sub create_oligo_annotation_files {
    my ($dbh, $fhLog, $sOutDir, $rAnnotMap, @anAnnotations) = @_;
    my @asAnnotationFiles;

    # Retrieve data for each annotation selected
    my $sSQL = join(' ', 'SELECT',
		    'r."reporterId" as "reporterId",', 
		    'rb."start" as "matchstart",',
		    'rb."end" as "matchend",',
		    'rb."score" as "score",',
		    'rb."crosshybRank" as "crosshybRank",', 
		    'bi."chromosomeId" as "chromosomeID",',
		    'bi."chromosomeName" as "chromosomeName"',
		    'FROM reporterbiosequencematch rb',
		    'LEFT JOIN biosequence bi ON rb."bioSequenceId"=bi.id',
		    'LEFT JOIN reporter r ON rb."reporterId"=r.id',
		    'WHERE rb."reporterLibraryAnnot"=?',
		    'ORDER BY r."id", rb."crosshybRank" DESC');
    my $oQuery = $dbh->prepare($sSQL);
    foreach my $nAnnotation (@anAnnotations){
	
	# Execute query
	my $n = $oQuery->execute($nAnnotation);
	die "Nothing found \n" if ( $n eq '0E0' );
	
	# Format output
	my %hColorCounter;
	my %hCrossHybOligos;
	my @aaOutput;
	$fhLog->print("\nPreparing oligo track '$rAnnotMap->{$nAnnotation}'\n");
	while (my @asRow = $oQuery->fetchrow_array() ){
	    my ($sID, $nStart, $nStop, $nScore, $nRank, $nChrID, $sChrName) = @asRow;
	    my $sChrID = join('', 'chr', $nChrID);
	    if ($nScore>=$ENV{OLIGO_SCORE_CUTOFF}){
		if ($nRank>$ENV{MAX_CROSSHYB_COUNT}){
		    print STDERR join("\t", $sID, $nScore, $nRank), "\n";
		}
	    
		$hCrossHybOligos{$sID}++ if ($nRank>$ENV{MAX_CROSSHYB_COUNT});
		push @aaOutput, [$sID, $sChrID, $rAnnotMap->{$nAnnotation}, 'oligo', $nStart, $nStop, $nScore, '.', '.', "$sID.$nRank"];
	    }
	}
	    
	# Print oligo's to file, make sure to exclude those oligos that have too many crosshybs
	my $sOutFilename = join('', $rAnnotMap->{$nAnnotation}, '.gff');
	$sOutFilename =~ s/ /-/g;
	open  OUTPUT, ">$sOutDir/$sOutFilename" or die "Can't open output file $sOutFilename: $!\n";
	foreach my $rRow (@aaOutput){
	    my @asRow = @$rRow;
	    my $sID   = shift @asRow;
	    print OUTPUT join("\t", @asRow), "\n" unless (exists($hCrossHybOligos{$sID}));
	}
	close OUTPUT;
	
	
	# Print some log info
	my $nDroppedOligos = scalar(keys(%hCrossHybOligos));
	$fhLog->print("   Features dropped because of crosshybs:    $nDroppedOligos\n");
	push @asAnnotationFiles, $sOutFilename;
    
    }
    $oQuery->finish();
    
    # Append retrieved
    return(@asAnnotationFiles);
}


# append_base_quickload_files (ContentString)
# 
# Creates the das_servers.txt, synonyms.txt and contents.txt files
sub append_base_quickload_files {
    my ($sOrganism, $fhLog, $sContentsString) = @_;
    my (@asFullname) = split /_/, $sContentsString;
    my $sOrgFull  = $ENV{ORGANISM}{$sOrganism}{name};
    my $sFullname = "$sOrgFull, Ensembl release $asFullname[2]-$asFullname[3], $asFullname[1]";
    
    `touch das_servers.txt`;
    `touch synonyms.txt`;
    open CONTENTS, ">>contents.txt" or die "Can't open output file 'contents.txt': $!\n";
    print CONTENTS join("\t", $sContentsString, $sFullname), "\n";
    close CONTENTS;
    $fhLog->print("\nAppended '$sFullname' to 'contents.txt' file");
}


# create_annots_file (ContentString)
# 
# Creates the annots.txt file
sub create_annots_file {
    my ($fhLog, $sOutdir, @asAnnotations) = @_;

    open ANNOTS, ">$sOutdir/annots.txt" or die "Can't open output file for writing 'annots.txt': $!\n";
    foreach my $sAnnotation (@asAnnotations){
	print ANNOTS $sAnnotation, "\n";
    }
    close ANNOTS;
    $fhLog->print("\nWrote annots.txt file\n");
}



# create_genome_seq_and_annot_files
#
# Subroutine that retrieves sequence data for all chromosomes and stores them as .bnib files
# also retrieves all gene features from the database
sub create_genome_seq_and_annot_files {
    my ($fhLog, $sOrganism, $sEnsBuild, $sOutdir) = @_; 
    my $sCore = $ENV{ORGANISM}{$sOrganism}{core};
    my @aoGeneObjects;

    # connect to ensembl database
    my $dbhEns = Bio::EnsEMBL::DBSQL::DBAdaptor->new( -user   => $ENV{ENSEMBL_USER},
						      -dbname => join('_', $sCore, $sEnsBuild),
						      -host   => $ENV{ENSEMBL_HOST},
						      -driver => "mysql",
						      -group  => "core");
    if (!$dbhEns){
	die "Failed to connect to database\n";
    }
    else{
	$fhLog->print("Connected to the ensembl database\n");
    }

    my $oSliceAdaptor = $dbhEns->get_SliceAdaptor;
    my $rSlices       = $oSliceAdaptor->fetch_all('chromosome');
    $fhLog->print("No chromosomes found in ensembl db!\n") if (scalar(@$rSlices) == 0);

    # fetch fasta sequences and store it in a separate file per chromosome
    #foreach my $oChr ( $rSlices->[0] ) {
    foreach my $oChr ( @$rSlices ) {
	my $sSeqname   = join('', 'chr', $ENV{ORGANISM}{$sOrganism}{chrId}{$oChr->seq_region_name()});
	my $nChrSize   = $oChr->end();
	my $sFastaname = join('', $sSeqname, '.fa');
	my $sBnibName  = join('', $sSeqname, '.bnib');

	# Add chromosome to mod_chromInfo.txt
	open CHR, ">>$sOutdir/mod_chromInfo.txt" or die "Can't open output file 'mod_chromInfo.txt': $!\n";
	print CHR join("\t", $sSeqname, $nChrSize), "\n";
	close CHR;
	
	# write sequences to file
	$fhLog->print("Writing sequence data for chromosome $sSeqname to file $sFastaname\n");
	my $sSequence = $oChr->seq();
	my $fhOut = new FileHandle;
	$fhOut->open(">$sOutdir/$sFastaname") or die "Can't open fasta file for writing: $!\n";
	$fhOut->print(">$sSeqname\n$sSequence\n");
	$fhOut->close();
	$fhLog->print("Appending sequence data for chromosome $sSeqname to file $ENV{COMPLETE_GENOME_SEQ}\n");
	my $fhGenome = new FileHandle;
	$fhGenome->open(">>$sOutdir/$ENV{COMPLETE_GENOME_SEQ}") or die "Can't open $ENV{COMPLETE_GENOME_SEQ} for writing: $!\n";
	$fhGenome->print(">$sSeqname\n$sSequence\n");
	$fhGenome->close();
	
	# convert fasta file to binary format
	`$ENV{BNIB_COMMAND} $sSeqname $sOutdir/$sFastaname $sOutdir/$sBnibName 1 >/dev/null`;
	`rm $sOutdir/$sFastaname`;
	
	# Get all genes on the chromosome
	my $rGenes = $oChr->get_all_Genes();
	$fhLog->print(join(' ', "Found", scalar(@$rGenes), "genes on chromosome $sSeqname\n\n"));
	open ANNOT, ">>$sOutdir/$ENV{GENE_ANNOTS_FILE}" or die "Can't write annotations to '$sOutdir/ensembl-genes.txt': $!\n";
	foreach my $oGene (@$rGenes){
	    #my $sID     = $oGene->external_name();
	    my $sID        = $oGene->stable_id();# unless ($sID);
	    my $rExons  = $oGene->get_all_Exons();
	    my $sStrand = $oGene->strand();
	    if ($sStrand == 1) {$sStrand = '+';}
	    else               {$sStrand = '-';}
	    
	    foreach my $oExon (@$rExons){
		my $nStart = $oExon->start();
		my $nEnd   = $oExon->end();
		print ANNOT join("\t", $sSeqname, 'orf_uncharacterized', 'exon', $nStart, $nEnd, '.', $sStrand, '.', $sID), "\n";
	    }
	}
	close ANNOT;

    }
    return ('ensembl-genes.gff');
}



# sdg_reannotate_and_peakcluster
#
# subroutine that retrieves current gene ontology files from www.geneontology.org,
# converts them to tab delimited files and stores goDivisions by goId in a hash
sub sdg_reannotate_and_peakcluster {
    my ($fhLog, $sOrganism, $sOutdir, $sGenomeAnnotationFile, $sGenespringMappings) = @_;

    # download sgd annotation files
    unless(-e "$ENV{TMP_DIR}") { mkdir($ENV{TMP_DIR}) or die "Could not create output dir $ENV{TMP_DIR}: $!\n"};

    my $sUrl     = $ENV{ORGANISM}{$sOrganism}{synonymsurl};
    my $sFile    = $ENV{ORGANISM}{$sOrganism}{synonymsfile};
    my $sOutFile = $ENV{TMP_DIR} . $sFile;

    # download sgd file
    $fhLog->print("\nDownloading SGD annotations file\n");
    if ( -f $sOutFile && ! -z $sOutFile ) {
#	return 1;
    } else {
	my $return = `wget -O $sOutFile $sUrl$sFile 2>&1`;
	my $ok = "\`$sOutFile\' saved";
	if ( ($return !~ /$ok/) ){ # error
	    die "Error saving SGD feature file $sOutFile:\n$return";
	}
    }
    
    # Parse sgd file
    my %hAnnotationLUT;
    $fhLog->print("Parsing SGD annotations file\n");
    open INPUT, "$sOutFile" or die "Can't open sgd feature file: $!";
    while (<INPUT>){
	chomp;
	my (@asLine) = split /\t/;
	if ($asLine[3]){
	    my $sSGDid             = uc($asLine[0]);
	    my $sSystematicName    = lc($asLine[3]);
	    my $sStandardName      = $asLine[4];
	    my $sClassification    = lc($asLine[1]);
	    my $sSubclassification = lc($asLine[2]);
	    $sClassification       = join('_', $sClassification, $sSubclassification) if ($sClassification eq 'orf');
	    $sStandardName         = '' if ($sClassification eq 'trna');
	    $hAnnotationLUT{$sSystematicName} = [$sSGDid, $sClassification, $sStandardName];
	}
    }
    close INPUT;
    
    # Read the genespring mappings file into a LUT
    my %hGenespringLUT;
    open GENESPRING, "$ENV{GENESPRING_MAPPINGS_DIR}/$sGenespringMappings" or die "Can't open genespring mappings file: $!\n";
    while (<GENESPRING>){
	s/\n//g; s/\r//g;
	next if (/^\s*$/);
	my ($sSGDid, $sOligoID)  = split /\t/;
	$hGenespringLUT{uc($sSGDid)} = $sOligoID;
    }
    close GENESPRING;
    
    
    # Reannotate the genome file and also write a peakcluster file
    open ANNOT, "$sOutdir/$sGenomeAnnotationFile" or die "Can't open annotation file: $!\n";
    open TMP, ">$sOutdir/$sGenomeAnnotationFile.tmp" or die "Can't open temp annotation file for writing: $!\n";
    open PEAKCLUSTER, ">$sOutdir/$ENV{PEAK_CLUSTER_FILE}" or die "Can't open peakcluster file for writing: $!\n";
    while(<ANNOT>){
	s/\n//g; s/\r//g;
	next if (/^\s*$/);
	my (@asAnnotation) = split /\t/;
	if (exists($hAnnotationLUT{lc($asAnnotation[8])}) ){
	    my ($sSGDid, $sClassification, $sStandardname) = @{$hAnnotationLUT{lc($asAnnotation[8])}};
	    $asAnnotation[8] = $sStandardname   if ($sStandardname);
	    $asAnnotation[1] = $sClassification if ($sClassification);
	    
	    # Write info to peakcluster file
	    if (exists($hGenespringLUT{$sSGDid})){
		print PEAKCLUSTER join("\t", $hGenespringLUT{$sSGDid}, $asAnnotation[1], $asAnnotation[0], $asAnnotation[3], $asAnnotation[4], $asAnnotation[6]), "\n";
	    }
	}
	print TMP join("\t", @asAnnotation), "\n";
    }
    close TMP;
    close PEAKCLUSTER;
    close ANNOT;
    `mv $sOutdir/$sGenomeAnnotationFile.tmp $sOutdir/$sGenomeAnnotationFile`;
}



# map_harbison_transfacs
#
# subroutine that maps putative transcription binding sites defined by the harbison lab
# and creates gff formatted mappings for each set
sub map_harbison_transfacs {
    my ($fhLog, $sOrganism, $sOutDir) = @_;
    my @asReturn;
    
    $fhLog->print("\nBlatting Harbison data to the new genome\n");
    my @asMappingFiles = @{$ENV{ORGANISM}{$sOrganism}{harbisonmappings}};
    foreach my $sMappingFile (@asMappingFiles){
	my %hLUT;
	my $sGffOut = basename($sMappingFile);
	`$ENV{BLAT_BIN} $sOutDir/$ENV{COMPLETE_GENOME_SEQ} $ENV{SUPPORTING_DATA_PATH}/$sMappingFile -out=blast8 $sOutDir/blat-output.txt`;
	
	my $fhBlat = new FileHandle;
	my $fhGff  = new FileHandle;
	$fhBlat->open("$sOutDir/blat-output.txt") or die "Couldn't open blat output file: $!\n";
	$fhGff->open(">$sOutDir/$sGffOut.gff") or die "Couldn't open mapping output file '$sGffOut.gff': $!\n";
	while (<$fhBlat>){
	    next if (/^\s*$/);
	    s/[\n\r]//g;
	    my ($sQdesc, $sHchr, $sScore, $nQLen, $sA, $sB, $nQstart, $nQstop, $nHstart, $nHstop, $sC, $sD) = split /\t/;
	    my ($sFactor, $sQchr, $nOstart, $nOstop, $nOffset, $nOlength, $sStrand) = split /\|/, $sQdesc;
	    if ( ($nOlength == $nQstop) and ($nQstart == 1)){
		$hLUT{$sFactor}++;
		my $sName  = join('.', $sFactor, $hLUT{$sFactor});
		my $nStart = $nHstart + $nOffset;
		my $nStop  = $nHstop  - $nOffset;
		$fhGff->print(join("\t", $sHchr, $sGffOut, 'transfac', $nStart, $nStop, '.', $sStrand, '.', $sName), "\n");
	    }
	}
	$fhBlat->close();
	$fhGff->close();
	push @asReturn, "$sGffOut.gff";
    }
    
    return (@asReturn);
}



# map_genome_wide_transcription
#
# subroutine that maps genome-wide transcription defined by David et al
# and creates sgr formatted mappings
sub map_genome_wide_transcription {
    my ($fhLog, $sOrganism, $sOutDir) = @_;
    my @asReturn;
    
    $fhLog->print("\nBlatting David data to the new genome\n");
    my @asMappingFiles = @{$ENV{ORGANISM}{$sOrganism}{transcriptionmappings}};
    foreach my $sMappingFile (@asMappingFiles){
	my %hLUT;
	my $sSgrOut = basename($sMappingFile);
	`$ENV{BLAT_BIN} $sOutDir/$ENV{COMPLETE_GENOME_SEQ} $ENV{SUPPORTING_DATA_PATH}/$sMappingFile -out=blast8 $sOutDir/blat-output.txt`;
	
	my $fhBlat = new FileHandle;
	my $fhSgrPos  = new FileHandle;
	my $fhSgrNeg  = new FileHandle;
	$fhBlat->open("$sOutDir/blat-output.txt") or die "Couldn't open blat output file: $!\n";
	$fhSgrPos->open(">$sOutDir/$sSgrOut-top.sgr") or die "Couldn't open mapping output file '$sSgrOut-top.sgr': $!\n";
	$fhSgrNeg->open(">$sOutDir/$sSgrOut-bot.sgr") or die "Couldn't open mapping output file '$sSgrOut-bot.sgr': $!\n";
	while (<$fhBlat>){
	    next if (/^\s*$/);
	    s/[\n\r]//g;
	    my ($sQdesc, $sHchr, $sScore, $nQLen, $sA, $sB, $nQstart, $nQstop, $nHstart, $nHstop, $sC, $sD) = split /\t/;
	    my ($sQchr, $nPos, $nValue, $sStrand, $nOlength) = split /\|/, $sQdesc;
	    if ( ($nOlength == $nQstop) and ($nQstart == 1)){
		my $nHpos = int( ($nHstart+$nHstop)/2 );
		if ($sStrand eq '+'){ $fhSgrPos->print(join("\t", $sHchr, $nHpos, $nValue), "\n"); }
		else                { $fhSgrNeg->print(join("\t", $sHchr, $nHpos, $nValue), "\n"); }
	    }
	}
	$fhBlat->close();
	$fhSgrPos->close();
	$fhSgrNeg->close();
	push @asReturn, "$sSgrOut.sgr";
    }
    
    return (@asReturn);
}




# retrieve_sgd_regulatory_mappings
#
# subroutine that downloads the gff file with transcription factor binding site positions
# from SGD and saves it as the file 'sgd_regulatory.gff'
sub retrieve_sgd_regulatory_mappings {
    my ($fhLog, $sOrganism, $sOutdir) = @_;
    my $sUrl = $ENV{ORGANISM}{$sOrganism}{sgdregulatory};
    my $sOut = join('/', $sOutdir, 'sgd_regulatory.tmp');
    $fhLog->print("Downloading transcription factor binding site mappings from SGD\n");
    
    my $return = `wget -O $sOut $sUrl 2>&1`;
    my $ok = "\`$sOut\' saved";
    if ( ($return !~ /$ok/) ){ # error
	die "Error saving SGD transcription factor binding site mappings file $sOut:\n$return";
    }
    
    # Open the file and replace roman numerals by integers
    my $fhTmpFile = new FileHandle;
    my $fhOutFile = new FileHandle;
    $fhTmpFile->open("$sOutdir/sgd_regulatory.tmp") or die "Can't open temporary SGD transfac mappings file: $!\n";
    $fhOutFile->open(">$sOutdir/sgd_regulatory.gff") or die "Can't open SGD transfac mappings file for writing: $!\n";
    while (<$fhTmpFile>){
	next if (/^#/);
	next if (/^\s*$/);
	my ($sChr, $sTrackName, @asRest) = split /\t/;
	next if ($sTrackName eq 'SGD');
	my $sOriginal = $sChr;
	$sChr       =~ s/chr//;
	$sTrackName = 'sgd_regulatory';
	if (exists($ENV{ORGANISM}{$sOrganism}{chrId}{$sChr})){
	    my $sNewChr = join('', 'chr', $ENV{ORGANISM}{$sOrganism}{chrId}{$sChr});
	    $fhOutFile->print(join("\t", $sNewChr, $sTrackName, @asRest));
	}
	else {
	    print STDERR "Warning: chromosome '$sOriginal' cannot be mapped onto IGB genome, skipping!\n";
	}
    }
    $fhTmpFile->close();
    $fhOutFile->close();
    `rm $sOutdir/sgd_regulatory.tmp`;
    return ('sgd_regulatory.gff');
}



# copy_lexA_control_plasmid_mappings
#
# subroutine that copies the lexAbirA plasmid mappings to the new genome
# the plasmid is always the same, so can be copied as a whole
sub copy_lexA_control_plasmid_mappings {
    my ($fhLog, $sOrganism, $sOutdir) = @_;
    my $sLexAbirA = $ENV{ORGANISM}{$sOrganism}{lexAcontrol};

    my $sGff  = "$sLexAbirA.gff";
    my $sBnib = "$sLexAbirA.bnib";
    $fhLog->print("\nAdding the LexA control plasmid mappings\n");
    copy("$ENV{SUPPORTING_DATA_PATH}/$sGff",  "$sOutdir/$sGff")  or die "Could not copy $ENV{SUPPORTING_DATA_PATH}/$sGff file to $sOutdir/$sGff: $!\n";
    copy("$ENV{SUPPORTING_DATA_PATH}/$sBnib", "$sOutdir/$sBnib") or die "Could not copy $ENV{SUPPORTING_DATA_PATH}/$sBnib bnib file to $sOutdir/$sBnib: $!\n";
    return ($sGff);
}



# map_rickyoung_signaling
#
# subroutine that creates the .sgr files for the
# Rick Young signaling data
sub map_rickyoung_signaling {
    my ($fhLog, $sOrganism, $sOutdirSgr, $sEnsDate, $reps, @anAnnotations) = @_;
    my $sRYSignaling = $ENV{SUPPORTING_DATA_PATH}."/RickYoung_SignalingData_UMCUids.txt";

    $fhLog->print("\nCreating .sgr files for Rick Young signaling data\n");

    foreach my $anAnnotation (@anAnnotations){

	# see how many data columns there are
	my @headerColumns = [];
	open (RYFILE, $sRYSignaling);
	while (<RYFILE>){
	    chomp $_;
	    @headerColumns = split '\t', $_;
	    last;
	}
	close RYFILE;
	
	# create the files
	my $hcnt = 0;
	while ( $hcnt < scalar(@headerColumns) ){
	    $hcnt++;
	    if ( $hcnt < 5 ){
		next;
	    }
	    
	    my $name = undef;
	    my $sgrfile = undef;
	    
	    open (RYFILE, $sRYSignaling);
	    my $lcnt = 0;
	    
	    while (<RYFILE>){
		chomp $_;
	    $lcnt++;
		
		if ( $_ =~ /^UMCU/ ){
		    if ( $hcnt > 5 ){
			close SGRFILE; # close previous file
		    }
		    $name = $headerColumns[($hcnt-1)];
		    $sgrfile = $sOutdirSgr . "/" . $name."_".$sEnsDate.".sgr"; 
		    $sgrfile =~ s/\s//g;
		    open (SGRFILE, ">$sgrfile") or die "unable to open sgr file $sgrfile for $name for writing\n";
		    next;
		}
		
		$_ =~ s/FLAG//g; # remove all 'FLAG's
		my @columns = split '\t', $_;
		my $chr = int($columns[2]);
		my $rep = $columns[0];
		
		my $repStart = $reps->{ $anAnnotation }{ $rep }{ matchstart };
		my $repEnd = $reps->{ $anAnnotation }{ $rep }{ matchend };
		if ( $repStart > $repEnd ){
		    $repStart = $reps->{ $anAnnotation }{ $rep }{ matchend };
		    $repEnd = $reps->{ $anAnnotation }{ $rep }{ matchstart };
		}

		my $newPos = int( ( ( $repEnd - $repStart ) / 2 ) + $repStart );

		if ( ! ( $newPos ) ){  # reporter is unmapped
		    next;
		}

		my $signal = $columns[($hcnt-1)];

		if ( ! $signal ){
		    $signal = 0; # IGB can't handle empty fields
		}

		print SGRFILE "chr" . $chr . "\t" . $newPos . "\t" . $signal . "\n"; 

	    } # while (<RYFILE
	    close RYFILE;
	} # while ( $hcnt < scalar(@headerColumns
    } # foreach my $nAnnotation
} # sub map_rickyoung_signaling



# retrieve_reporter_startstop
#
# retrieves reporter start and stop positions and stores it in memory
sub retrieve_reporters {
    my ($dbh, $fhLog, $sOrganism, @anAnnotations) = @_;
    my $reps = {};

    my $sSQL = join(' ', 'SELECT',
		    'r."reporterId" as "reporterId",', 
		    'rbm."start" as "matchstart",',
		    'rbm."end" as "matchend",',
		    'rbm."score" as "score",',
		    'rbm."crosshybRank" as "crosshybRank"', 
		    'FROM reporter r, reporterbiosequencematch rbm',
		    'WHERE rbm."reporterId" = r.id',
		    'AND ( rbm."crosshybRank" = 1 OR rbm."crosshybRank" IS NULL )',
		    'AND "reporterLibraryAnnot" =?');
    my $oQuery = $dbh->prepare($sSQL);
    foreach my $nAnnotation (@anAnnotations){
	
	# Execute query
	my $n = $oQuery->execute($nAnnotation);
	die "Nothing found \n" if ( $n eq '0E0' );

	while (my @asRow = $oQuery->fetchrow_array() ){
	    my ($sID, $nStart, $nStop, $nScore, $nRank) = @asRow;
	    $reps->{ $nAnnotation }{ $sID } = { matchstart => $nStart,
				matchend => $nStop,
				score => $nScore,
				crosshybrank => $nRank,
			    };
	}
    }

    return $reps;

} # sub retrieve_reporters
