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
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $List, $Trusted, $MainPath, $MolType, $OutPath);

$Usage = "\tUsage: MakeBlastDBs.pl <Project_Name> <List_File.ext> <Trusted_ORFeome.fasta> <Main_Path> <Molecule type>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName    = $ARGV[0];
$List           = $ARGV[1];
$Trusted        = $ARGV[2];
$MolType        = $ARGV[3];  # nucl or prot
$OutPath        = $ARGV[4];

my ($Project, $ORFeomesPath, $BlastPath, $TrustedORFeomeDb, $SeqExt,
	$i, $n, $Qry, $InputFile, $Db, $cmd, $LogFile, $TrustedPrefix);
my (@List);

#$Project              = $OutPath ."/". $ProjectName;
$Project = $OutPath;
$BlastPath            = $OutPath ."/". "Blast";
$ORFeomesPath         = $OutPath ."/". "ORFeomes";

$TrustedPrefix        = Prefix($Trusted);
$Trusted              = $ORFeomesPath ."/". $Trusted;
$TrustedORFeomeDb     = $BlastPath ."/". $TrustedPrefix . "Db";
$LogFile              = $OutPath ."/". $ProjectName . ".log";

if ($MolType eq "nucl"){
	$SeqExt = ".ffn";
}elsif($MolType eq "prot"){
	$SeqExt = ".faa";
}

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

print "\nBuilding data bases:\n";

MakeDir($BlastPath); 

@List = ReadFile($List); 
$n = scalar@List;

$cmd = `makeblastdb -in $Trusted -dbtype $MolType -parse_seqids -out $TrustedORFeomeDb`;

for ($i=0; $i<$n; $i++){
	$Qry = $List[$i];
	$InputFile = $ORFeomesPath ."/". $Qry . $SeqExt;     
	$Db = $BlastPath ."/". $Qry;
	$cmd = `makeblastdb -in $InputFile -dbtype $MolType -parse_seqids -out $Db`;  
	Progress($n, $i);
}
exit;