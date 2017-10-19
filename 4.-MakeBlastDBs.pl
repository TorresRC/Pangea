#!/usr/bin/perl -w
#################################################################################
#Scipt ContigsLength.pl                                                         #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          10 de abril de 2017                                             #
#################################################################################
use strict;                                                                      
use lib '/Users/rc/lib';                              
use Routines;                                                                    

my ($Usage, $ProjectName, $List, $TrustedORFeome);

$Usage = "\tUsage: MakeBlastDBs.pl <Project_Name> <List_File.ext> <Trusted_ORFeome.fasta>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List = $ARGV[1];
$TrustedORFeome = $ARGV[2];

my ($MainPath, $Project, $ORFeomesPath, $MainList, $BlastPath, $TrustedORFeomeDb, $SeqExt,
	$i, $n, $Qry, $InputFile, $Db, $cmd, $LogFile, $TrustedORFeomePrefix);
my (@List);                                                     

$MainPath             = "/Users/rc/CoreGenome";
$Project              = $MainPath ."/". $ProjectName;
$MainList             = $Project ."/". $List;
$BlastPath            = $MainPath ."/". "Blast";
$ORFeomesPath         = $MainPath ."/". "ORFeomes" ."/". "Sorted" ."/". "Filtered";
$TrustedORFeomePrefix = Prefix($TrustedORFeome);
$TrustedORFeome       = $ORFeomesPath ."/". $TrustedORFeome;
$TrustedORFeomeDb     = $BlastPath ."/". $TrustedORFeomePrefix . "Db";
$SeqExt               = ".ffn";
$LogFile              = $Project ."/". $ProjectName . ".log";

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

print "\nBuilding data bases:\n";

MakeDir($BlastPath);                                                             

@List = ReadFile($MainList);                                                     
$n = scalar@List;

$cmd = `makeblastdb -in $TrustedORFeome -dbtype nucl -parse_seqids -out $TrustedORFeomeDb`;

for ($i=0; $i<$n; $i++){                                                         
	$Qry = $List[$i];                                                        

	$InputFile = $ORFeomesPath ."/". $Qry . $SeqExt;                             
	$Db = $BlastPath ."/". $Qry;                                          

	$cmd = `makeblastdb -in $InputFile -dbtype nucl -parse_seqids -out $Db`;  

	Progress($n, $i);
}
exit;
