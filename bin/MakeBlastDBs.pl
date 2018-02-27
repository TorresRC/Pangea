#!/usr/bin/perl -w
#################################################################################
#Scipt MakeBlastDBs.pl                                                          #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          10 de abril de 2017                                             #
#################################################################################
use strict;
use FindBin;
use lib "$FindBin::Bin/../lib_nt";
use Routines;

my ($Usage, $ProjectName, $List, $TrustedORFeome, $MainPath, $MolType);

$Usage = "\tUsage: MakeBlastDBs.pl <Project_Name> <List_File.ext> <Trusted_ORFeome.fasta> <Main_Path> <Molecule type>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$MainPath       = $ARGV[0];
$ProjectName    = $ARGV[1];
$List           = $ARGV[2];
$TrustedORFeome = $ARGV[3];
$MolType        = $ARGV[4];  # nucl or prot

my ($Project, $ORFeomesPath, $MainList, $BlastPath, $TrustedORFeomeDb, $SeqExt,
	$i, $n, $Qry, $InputFile, $Db, $cmd, $LogFile, $TrustedORFeomePrefix);
my (@List);

$Project              = $MainPath ."/". $ProjectName;
$MainList             = $Project ."/". $List;
$BlastPath            = $Project ."/". "Blast";
$ORFeomesPath         = $Project ."/". "ORFeomes";
$TrustedORFeomePrefix = Prefix($TrustedORFeome);
$TrustedORFeome       = $ORFeomesPath ."/". $TrustedORFeome;
$TrustedORFeomeDb     = $BlastPath ."/". $TrustedORFeomePrefix . "Db";

$LogFile              = $Project ."/". $ProjectName . ".log";

if ($MolType eq "nucl"){
	$SeqExt = ".ffn";
}elsif($MolType eq "prot"){
	$SeqExt = ".faa";
}

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

print "\nBuilding data bases:\n";

MakeDir($BlastPath); 

@List = ReadFile($MainList); 
$n = scalar@List;

$cmd = `makeblastdb -in $TrustedORFeome -dbtype $MolType -parse_seqids -out $TrustedORFeomeDb`;

for ($i=0; $i<$n; $i++){
	$Qry = $List[$i];

	$InputFile = $ORFeomesPath ."/". $Qry . $SeqExt;     
	$Db = $BlastPath ."/". $Qry;

	$cmd = `makeblastdb -in $InputFile -dbtype $MolType -parse_seqids -out $Db`;  

	Progress($n, $i);
}
exit;
