#!/usr/bin/perl -wT


# Cgi script used to analyze the path_blast database
# All HTML generated in script file


# MODULES

use strict;
use FileHandle;
use CGI;



#----------------------------------------------------------------------------
# start_igb_html($sTitle, $sJscript)
#
# Print the page header of the shops output page with correct scripting stuff
sub start_igb_html{
    my $sTitle   = shift;
    my $sHeader = '';

    $sHeader  = &header(-type=>'text/html', -expires=>'now');
    $sHeader .= &start_html(-title   =>$sTitle,
			    -author  =>$ENV{MAILTO},
			    -meta    =>{'keywords'=>'SHOPS'},
			    -style   =>{'src'=>$ENV{STYLE_URL}},
			    -head    =>join('',
					    "<script language='JavaScript' src='$ENV{JSCRIPT}overlib.js' type='text/javascript'></script>",
					    "<script language='JavaScript' src='$ENV{JSCRIPT}overlib_shadow.js' type='text/javascript'></script>",
					    "<script language='JavaScript' src='$ENV{JSCRIPT}overlib_exclusive.js' type='text/javascript'></script>",
					    "<script language='JavaScript' src='$ENV{JSCRIPT}overlib_centerpopup.js' type='text/javascript'></script>",
					    "<script language='JavaScript' src='$ENV{JSCRIPT}overlib_followscroll.js' type='text/javascript'></script>",
					    "<script language='JavaScript' src='$ENV{JSCRIPT}overlib_hideform.js' type='text/javascript'></script>",
					    '<link rel="icon" href="/tools/shops/favicon.ico" type="image/x-icon">',
					    '<link rel="shortcut icon" href="/tools/shops/favicon.ico" type="image/x-icon">'),
			    -BGCOLOR =>'white',
			    -onunload=>'nd();');
    $sHeader .= "<div id='overDiv' style='position:absolute; visibility:hidden; z-index:1000;'></div>\n";
    $sHeader .= "<H1>$sTitle</H1>\n";
    
    return ($sHeader);
}



#----------------------------------------------------------------------------
# help_popup_hash
#
# Returns javascript code for a popup window with a help message
sub help_popup_hash{
    my $sKey     = shift;
    my $sMessage = '- empty -';
    my $nWidth   = $ENV{HELP_MESSAGE_WIDTH};
    my %hReturn;

    $sMessage = $ENV{HELP_MESSAGES}{$sKey} if (exists($ENV{HELP_MESSAGES}{$sKey}));
    $hReturn{'href'}        = '#';
    $hReturn{'onMouseOver'} = "overlib('$sMessage', WIDTH, $nWidth, BORDER, 1, BGCOLOR, '$ENV{HELP_MESSAGE_BGCOLOR}', CAPCOLOR, '#000000', FGCOLOR, '$ENV{HELP_MESSAGE_FGCOLOR}');";
    $hReturn{'onMouseOut'}  = 'nd();';
    return (\%hReturn);
}



#----------------------------------------------------------------------------
# scrolling_list_genomes
#
# Provides a scrolling list element with all species in it
# as a second parameter it returns a list with non-redundant species!
sub scrolling_list_genomes{
    my @asValues;
    my %hsLabels;
    my $sContents = join('/', $ENV{IGB_QUICKLOAD_PATH}, 'contents.txt');
    
    # Test if quickload file exists
    &errorMessage("Could not find a quickload directory at the specified location", $ENV{STYLE_URL}, $ENV{MAILTO})
	unless (-f $sContents);
	
    # Open the quickload contents.txt file
    open CONTENTS, "<$sContents" or &errorMessage("Could not open contents.txt file!", $ENV{STYLE_URL}, $ENV{MAILTO});
    while(<CONTENTS>){
	next if (/^\s*$/);
	s/\n//g; s/\r//g;
	my ($sGenomeDir, $sGenomeDesc) = split /\t/;
	push @asValues, $sGenomeDir;
	$hsLabels{$sGenomeDir} = $sGenomeDesc;
    }
    close CONTENTS;

    # Generate the code for the box
    my $sList = &scrolling_list(-name     => 'QUICKLOAD_GENOME',
				-values   => [@asValues],
				-size     => 5,
				-multiple => 'false',
				-defaults => [$asValues[0]],
				-labels   => \%hsLabels);
				
    return ($sList);
}



#----------------------------------------------------------------------------
# get_oligo_positions($sGenomeDir)
#
# Get mappings for all oligo's
sub get_oligo_positions {
    my $sGenomeDir = shift @_;
    my %hOligoMappings;

    # Get names of oligo tracks
    my @asTracks;
    my $sAnnotsFile = join('/', $ENV{IGB_QUICKLOAD_PATH}, $sGenomeDir, 'annots.txt');
    open ANNOTS, "<$sAnnotsFile" or &errorMessage("Could not open annots.txt file!", $ENV{STYLE_URL}, $ENV{MAILTO});
    while(<ANNOTS>){
	next if (/^\s*$/);
	s/\n//g; s/\r//g;
	push @asTracks, $_ if (/features\.gff/);
    }
    close ANNOTS;
    
    # Read each oligo track and extract mappings (based on extensions!)
    foreach my $sTrack (@asTracks){
	my $sTrackFormat;
	my $sTrackFile = join('/', $ENV{IGB_QUICKLOAD_PATH}, $sGenomeDir, $sTrack);
	
	# Check format of the track file
	my $sCheckTrack = lc($sTrack);
	if ($sCheckTrack=~ /gff$/)    { $sTrackFormat = 'gff'; }
	elsif ($sCheckTrack =~ /bed$/){ $sTrackFormat = 'bed'; }
	elsif ($sCheckTrack =~ /sgr$/){ next; }
	else{ &errorMessage("Unknown file format for track file: $sTrackFile!", $ENV{STYLE_URL}, $ENV{MAILTO}); }
	
	# Read oligo mappings
	open TRACK, "<$sTrackFile" or &errorMessage("Could not open $sTrackFile!", $ENV{STYLE_URL}, $ENV{MAILTO});
	while(<TRACK>){
	    next if (/^\s*$/);
	    s/\n//g; s/\r//g;
	    s/\"//g; s/\'//g; # remove quotes!
	    my (@asLine) = split /\t/;
	    if ($sTrackFormat eq 'gff'){
		if (scalar(@asLine)>8){
		    my ($sOligo, $nVersion) = split /:/, uc($asLine[8]);
		    $hOligoMappings{$sOligo} = {} unless (exists($hOligoMappings{$sOligo}));
		    push @{$hOligoMappings{$sOligo}{$asLine[0]}}, [($asLine[3], $asLine[4])];
		}
	    }
	    if ($sTrackFormat eq 'bed'){
		if (scalar(@asLine)>3){
		    my ($sOligo, $nVersion) = split /:/, uc($asLine[3]);
		    $hOligoMappings{$sOligo} = {} unless (exists($hOligoMappings{$sOligo}));
		    push @{$hOligoMappings{$sOligo}{$asLine[0]}}, [($asLine[3], $asLine[4])];
		}
	    }
	}
    }
    
    return(\%hOligoMappings);
}


return 1;
