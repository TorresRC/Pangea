#!/usr/bin/perl -w
#################################################################################
#Scipt ContigsLength.pl                                                         #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          11 de octubre de 2017                                           #
#################################################################################
use strict; 
use FindBin;
use lib "$FindBin::Bin/lib";
use Routines;

my ($Usage, $ProjectName, $List, $TrustedORFeome, $MainPath);

$Usage = "\tUsage: FilterOrfeomes.pl <Project_Name> <List_File.ext> <Trusted_ORFeome.fasta> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List = $ARGV[1];
$TrustedORFeome = $ARGV[2];
$MainPath = $ARGV[3];

my ($Project, $ORFeomesPath, $MainList, $SeqExt, $Qry, $InputFile, $cmd, $LogFile, $FilteredORFsPath,
    $FilteredTrustedORFeome, $OutputFile, $ClstrsPath, $Clstr, $TrustedFile);
my ($i, $n);
my (@List);                                                     

$SeqExt = ".ffn";
$Project = $MainPath ."/". $ProjectName;
$MainList = $Project ."/". $List;
$ORFeomesPath = $MainPath ."/". "ORFeomes/Sorted";
$FilteredORFsPath = $ORFeomesPath ."/". "Filtered";
$ClstrsPath = $MainPath ."/". "GenesClusters";
$TrustedFile = $ORFeomesPath ."/". $TrustedORFeome;
$FilteredTrustedORFeome = $FilteredORFsPath ."/". $TrustedORFeome;

$LogFile = $Project ."/". $ProjectName . ".log";

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

MakeDir($FilteredORFsPath);
MakeDir($ClstrsPath);

@List = ReadFile($MainList);                                                     
$n = scalar@List;

print "\nFiltering Duplicates:\n";
$cmd = `cd-hit-dup -i $TrustedFile -o $FilteredTrustedORFeome -m false -e 0.20`;

for ($i=0; $i<$n; $i++){                                                         
	$Qry = $List[$i];

	$InputFile = $ORFeomesPath ."/". $Qry . $SeqExt;
        $OutputFile = $FilteredORFsPath ."/". $Qry . $SeqExt;
        $Clstr = $FilteredORFsPath ."/". "*.clstr";

	$cmd = `cd-hit-dup -i $InputFile -o $OutputFile -m false -e 0.05`;
        $cmd = `mv $Clstr $ClstrsPath`;
        
        Progress($n,$i);

} 
exit;
