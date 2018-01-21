#!/usr/bin/perl -w
#################################################################################
#Scipt ContigsLength.pl                                                         #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          10 de abril de 2017                                             #
#################################################################################
use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $List, $TrustedORFeome, $MainPath);

$Usage = "\tUsage: MakeBlastDBs.pl <Project_Name> <List_File.ext> <Trusted_ORFeome.fasta> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName    = $ARGV[0];
$List           = $ARGV[1];
$TrustedORFeome = $ARGV[2];
$MainPath       = $ARGV[3];

my ($Project, $ORFeomesPath, $MainList, $BlastPath, $TrustedORFeomeDb, $SeqExt,
	$i, $n, $Qry, $InputFile, $Db, $cmd, $LogFile, $TrustedORFeomePrefix);
my (@List);

$Project              = $MainPath ."/". $ProjectName;
$MainList             = $Project ."/". $List;
$BlastPath            = $Project ."/". "Blast";
$ORFeomesPath         = $Project ."/". "ORFeomes" ."/". "Sorted";
$TrustedORFeomePrefix = Prefix($TrustedORFeome);
$TrustedORFeome       = $ORFeomesPath ."/". $TrustedORFeome;
$TrustedORFeomeDb     = $BlastPath ."/". $TrustedORFeomePrefix . "Db";
$SeqExt               = ".faa";

$LogFile              = $Project ."/". $ProjectName . ".log";

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

print "\nBuilding data bases:\n";

MakeDir($BlastPath); 

@List = ReadFile($MainList); 
$n = scalar@List;

$cmd = `makeblastdb -in $TrustedORFeome -dbtype prot -parse_seqids -out $TrustedORFeomeDb`;

for ($i=0; $i<$n; $i++){
	$Qry = $List[$i];

	$InputFile = $ORFeomesPath ."/". $Qry . $SeqExt;     
	$Db = $BlastPath ."/". $Qry;

	$cmd = `makeblastdb -in $InputFile -dbtype prot -parse_seqids -out $Db`;  

	Progress($n, $i);
}
exit;
