#!/usr/bin/env perl

# Normalize affy cel files with tilecore.exe and wine

# MODULES
use strict;
use Getopt::Long;
use Cwd;

# GLOBALS
$ENV{WINE}       ||= "wine";
$ENV{TILECORE}   ||= "/home/hugheslab1/projects/bin/tilecore.exe";

# GET ARGUMENTS
my $sHelp         = 0;
my $sInput        = '';
my $sCelPath      = './';
my $sCelCtrlPath  = '';
my $sBpmapPath    = './';
my $sOutputSuffix = '';
my $sNormalize    = 1;
my $sExport       = 0;
my $nBandwidth    = 20;
my $nPvalScale    = 1;
my $nSigScale     = 2;
my $nTailType     = 2;
my $nNscale       = 1;
my $nTgtInt       = 100;
my $nOutputType   = 1;
my $nFormatType   = 1;
my $nPmonly       = 0;
my $nIndividual   = 0;
my $nNormTogether = 1;
GetOptions("help!"           => \$sHelp,
           "input:s"         => \$sInput,
           "cel_path:s"      => \$sCelPath,
           "cel_ctrl_path:s" => \$sCelCtrlPath,
           "bpmap_path:s"    => \$sBpmapPath,
           "normalize:i"     => \$sNormalize,
           "export:i"        => \$sExport,
           "band:i"          => \$nBandwidth,
           "pval_scale:i"    => \$nPvalScale,
           "sig_scale:i"     => \$nSigScale,
           "tail_type:i"     => \$nTailType,
           "nscale:i"        => \$nNscale, 
           "tgt:i"           => \$nTgtInt,
           "output_type:i"   => \$nOutputType,
           "output_suffix:s" => \$sOutputSuffix,
           "format_type:i"   => \$nFormatType,
           "pmonly!"         => \$nPmonly,
           "individual!"     => \$nIndividual,
           "nt:i"            => \$nNormTogether);


# PRINT HELP
$sHelp = 1 unless($sInput and $sCelPath and $sBpmapPath);
if ($sHelp) {
    die <<HELP

    tas_normalize.pl -input -cel_path -bpmap_path [optional parameters]

    Arguments:
    -input
      Tab-delimited input file with cel files to analyze
      row format: bpmap   control_name   ctrlCEL1;ctrlCEL2;   treatment_name   treatCEL1;treatCEL2;
      Alternatively, if you just want to normalize single CEL file(s), you can skip the last two
      columns:    bmap    cel_name       cel1;cel2;cel3
      Multiple CEL files can be normalized together in a group if they are separated by ';'.
    -cel_path
      Path where cel files can be found
    -bpmap_path 
      Path where the bpmap files can be found
      
    Optional:
    -cel_ctrl_path
      Path to where the control cel files can be found, if different from treatment
    -normalize
      normalization_flag (0 = false, 1 = true). Default: 1
    -export 
      export_as_text_flag (0 = false, 1 = true). Default: 0
    -band 
      bandwidth in nucleotides for averaging over probes. Default: 20
    -pval_scale 
      pvalue_scale (0 = linear, 1 = -10log10). Default: 1
    -sig_scale 
      signal_scale (0 = linear, 1 = log10, 2 = log2). Default: 2
    -tail_type 
      tail_type (0 = one sided lower, 1 = one sided upper, 2 = two sided). Default: 2
    -nscale
      boolean_scale_to_target_intensity (0 = false, 1 = true). Default: 1
    -tgt 
      normalization_target_intensity. Default: 100
      default: 100
    -nt
      Normalize together. Default: $nNormTogether
    -output_type 
      0=both signal and p-values, 1=signal only, 2=p-values only). Default: 1
    -output_suffix
      Optional suffix that is added to the output files. Default: none
    -format_type 
      0=CHP file, 1=BAR file). Default: 1
    -individual
      Add this flag to normalize individual cel files
    -pmonly
      Only use PM probes

HELP
}


###########
## START ##
###########

# Some simple conversions
$sCelPath     =~ s/\/$//;
$sCelCtrlPath =~ s/\/$//;
$sBpmapPath   =~ s/\/$//;

# Create syntax for common parameters
my $sCommonParameters = join(' ', "-n $sNormalize -export $sExport -band $nBandwidth -pval_scale $nPvalScale",
                                  "-sig_scale $nSigScale -tail_type $nTailType -nscale $nNscale",
                                  "-tgt $nTgtInt -output_type $nOutputType -format_type $nFormatType -nt $nNormTogether");
$sCommonParameters   .= " -pmonly" if ($nPmonly);


# Read the input file and store each job command in a hash
# The hash is used to avoid identical jobs from being submitted multiple times
my %hJobs;
open INPUT, $sInput or die "Can't open input file: $!\n";
while (<INPUT>){
   next if (/^\s*$/);
   next if (/^\s*#/);
   s/[\n\r]//g;
   my ($sBpmap, $sControlName, $sControlCel, $sTreatmentName, $sTreatmentCel) = split /\t/;
   my @asControlCel   = split /;/, $sControlCel;
   my @asTreatmentCel = split /;/, $sTreatmentCel;
   
   # Make sure that at least a name and control cel file are specified
   unless($sControlName and $sControlCel){
      print STDERR "Warning: Skipping line, no cel files defined\n";
      next;
   }
   
   # Add path to bpmap file
   $sBpmap = join('/', $sBpmapPath, $sBpmap);
   
   
   # See if there is a treatment, otherwise normalize cel file(s) individually 
   if (@asTreatmentCel){
      my $sOutput    = "${sTreatmentName}_vs_${sControlName}";
      $sOutput       = join('_', $sOutput, $sOutputSuffix) if ($sOutputSuffix);
      my $sControl   = $sCelCtrlPath ? &get_cel_syntax('control', $sCelCtrlPath, @asControlCel) : &get_cel_syntax('control', $sCelPath, @asControlCel);
      my $sTreatment = &get_cel_syntax('treatment', $sCelPath, @asTreatmentCel);
      my $sCommand   = join(' ', $ENV{WINE}, $ENV{TILECORE}, "-analysis", "-bpmap $sBpmap",
                                 $sControl, $sTreatment, "-out $sOutput", 
                                 "-type 1", $sCommonParameters);
      $hJobs{$sCommand}++;
                                   
      if ($nIndividual){
         $sOutput   = "$sControlName-single";
         $sOutput   = join('_', $sOutput, $sOutputSuffix) if ($sOutputSuffix);
         my $sCel   = $sCelCtrlPath ? &get_cel_syntax('cel', $sCelCtrlPath, @asControlCel) : &get_cel_syntax('cel', $sCelPath, @asControlCel);
         $sCommand  = join(' ', $ENV{WINE}, $ENV{TILECORE}, "-analysis", "-bpmap $sBpmap",
                                "-type 0", $sCel, "-out $sOutput", $sCommonParameters);
         $hJobs{$sCommand}++;

         $sOutput   = "$sTreatmentName-single";
         $sOutput   = join('_', $sOutput, $sOutputSuffix) if ($sOutputSuffix);
         $sCel      = &get_cel_syntax('cel', $sCelPath, @asTreatmentCel);
         $sCommand  = join(' ', $ENV{WINE}, $ENV{TILECORE}, "-analysis", "-bpmap $sBpmap",
                                "-type 0", $sCel, "-out $sOutput", $sCommonParameters);
         $hJobs{$sCommand}++;
      }
   }
   else{
      my $sOutput   = "$sControlName-single";
      $sOutput      = join('_', $sOutput, $sOutputSuffix) if ($sOutputSuffix);
      my $sCel      = &get_cel_syntax('cel',$sCelPath,  @asControlCel);
      my $sCommand  = join(' ', $ENV{WINE}, $ENV{TILECORE}, "-analysis", "-bpmap $sBpmap",
                                "-type 0", $sCel, "-out $sOutput", $sCommonParameters);
      $hJobs{$sCommand}++;
   }
   
   
}

# Now submit all jobs
foreach my $sJob (keys(%hJobs)){
   &submitjob($sJob);
   #print $sJob, "\n";
}


#################
## SUBROUTINES ##
#################

# get_cel_syntax($cel-type, @asCelfiles)
#
# Creates the cel file syntax for the tilecore commandline
sub get_cel_syntax {
   my ($sType, $sPath, @asCelFiles) = @_;
   my $sReturn;
   foreach my $sCelFile (@asCelFiles){
      $sReturn .= "-$sType $sPath/$sCelFile " if ($sCelFile);
   }
   return $sReturn;
}



# submitjob($sJobCommand)
#
# Submits the job to the node that has wine installed
sub submitjob {
   my ($sJobCmd) = @_;
   
   # Create pbs output dir if it doesn't exist yet
   my $sPbsOutDir = join ('/', $ENV{'HOME'}, 'pbs-output');
   unless (-e $sPbsOutDir) { mkdir($sPbsOutDir) or die "Could not create '$sPbsOutDir' for pbs job output: $!\n"};

   # Get current working directory
   my $sCwd     = getcwd;
   my $sJobPath = join(':', $ENV{'PATH'}, '/opt/sw/bin');
   
   open QSUB, "|qsub" or die "Can't execute qsub command: $!\n";
   print QSUB <<SUBMIT;
#PBS -S /bin/bash
#PBS -e localhost:$sPbsOutDir
#PBS -o localhost:$sPbsOutDir
#PBS -j oe
#PBS -l walltime=01:00:00,nodes=1:ppn=1,mem=1024M
#PBS -r n
#PBS -V
#PBS -N tas-normalize
cd $sCwd
export PATH='$sJobPath'
echo Torque run command: $sJobCmd
$sJobCmd
SUBMIT
   close QSUB;

}
