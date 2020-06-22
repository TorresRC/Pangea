#!/usr/bin/perl -w
#################################################################################
#Scipt MakeBlastDBs.pl                                                          #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          April 10th 2017                                                 #
#################################################################################
use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $Query, $MolType, $OutPath);

$Usage = "\tUsage: MakeBlastDBs.pl <Query name> <Molecule type> <OutPath>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$Query   = $ARGV[0];
$MolType = $ARGV[1];  # nucl or prot
$OutPath = $ARGV[2];

my ($BlastPath, $BlastDbPath, $SamplesFeaturesPath, $SeqExt, $InputFile, $Db,
	$cmd, $LogFile, $TrustedPrefix);

$LogFile = $OutPath ."/". "Errors.log";

$SamplesFeaturesPath = $OutPath   ."/". "DataFiles" ."/". "SamplesFeatures";
$BlastPath           = $OutPath   ."/". "DataFiles" ."/". "Blast";
$BlastDbPath         = $BlastPath ."/". "DBs";


if ($MolType eq "nucl"){
	$SeqExt = ".ffn";
}elsif($MolType eq "prot"){
	$SeqExt = ".faa";
}

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

MakeDir($BlastPath);
MakeDir($BlastDbPath); 

$InputFile = $SamplesFeaturesPath ."/". $Query . $SeqExt;     
$Db = $BlastDbPath ."/". $Query . "_BlastDb";
system("makeblastdb -in $InputFile -dbtype $MolType -parse_seqids -out $Db >/dev/null 2>&1");

close STDERR;
exit;