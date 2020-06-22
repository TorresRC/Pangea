#################################################################################
#Scipt FormatFeatures.pl                                                        #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          October 19 2017                                                 #
#################################################################################
use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use List::MoreUtils qw{any};
use Routines;

my ($Usage, $Query, $AnnotationPath, $MolType, $OutPath);

$Usage = "\tUsage: FormatFeatures.pl <Query List> <Annotation Path> <Molecule Type> <Outpath>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;

$Query          = $ARGV[0];
$AnnotationPath = $ARGV[1];
$MolType        = $ARGV[2];    #nucl or prot
$OutPath        = $ARGV[3];

my ($LogFile, $DataFiles, $SamplesFeaturesPath, $Ext, $Prefix, $OriginalFile,
	 $FormatedFile, $Line, $FormatedHeader);
my (@FeatureCollection, @Header);

$LogFile   = $OutPath ."/". "Errors.log";
$DataFiles = $OutPath ."/". "DataFiles";
$SamplesFeaturesPath = $DataFiles ."/". "SamplesFeatures";

if ($MolType eq "nucl"){
	$Ext = ".ffn";
}elsif($MolType eq "prot"){
	$Ext = ".faa";
}

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

MakeDir($SamplesFeaturesPath);

$Prefix = Prefix($Query);
$OriginalFile = $AnnotationPath ."/". $Prefix ."/". $Prefix.$Ext;
	
@FeatureCollection = ReadFile($OriginalFile);
$FormatedFile = $SamplesFeaturesPath ."/". $Prefix.$Ext;	
open(FILE, ">$FormatedFile");
foreach $Line (@FeatureCollection){	
	if ($Line =~ "^>" ){
		@Header = split (" ", $Line);
		$FormatedHeader = $Header[0];
		chomp $FormatedHeader;
		print FILE "$FormatedHeader\n";
	}else{
		$Line =~ s/-//g;
		print FILE "$Line\n";
	}
}
close FILE;

close STDERR;
exit;