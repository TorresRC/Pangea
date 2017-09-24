#!/usr/bin/perl -w                                                               

#################################################################################
#   Programa Make Blast DB                                                      #
#   Nota: Se debe ajustar la ruta de lib y de la variable $PathSeq a las que    #
#   realmente se tengan donde se instalción el programa.                        #
#                                                                               #
# Programador:   Roberto C. Torres                                              #
# Fecha:         10 de abril de 2017                                            #
#################################################################################
use strict;                                                                      
use lib '/home/bioinformatica/CoreGenome/src/lib';
#use lib '/Users/Roberto/Documents/lib';                              
use Routines;                                                                    

my ($Usage, $ProjectName, $List);

$Usage = "\tUsage: MakeBlastDBs.pl <Project Name> <List File Name>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List = $ARGV[1];                                                               

my ($MainPath, $Project, $ORFeomesPath, $MainList, $BlastPath, $ext,
	$i, $n, $Qry, $InputFile, $Db, $cmd, $LogFile);
my (@List);                                                     

#$MainPath = "/Users/Roberto/CoreGenome";
$MainPath = "/home/bioinformatica/CoreGenome";
$Project = $MainPath ."/". $ProjectName;
$MainList = $Project ."/". $List;                                               
#$ORFeomesPath = $Project ."/". "ORFeomes";
$ORFeomesPath = $MainPath ."/". "ORFeomes"; 
#$BlastPath = $Project ."/". "Blast";
$BlastPath = $MainPath ."/". "Blast";
$ext = ".fasta";
$LogFile = $Project ."/". $ProjectName . ".log";


open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

MakeDir($BlastPath);                                                             

@List = ReadFile($MainList);                                                     
$n = scalar@List;                                                                

for ($i=0; $i<$n; $i++){                                                         
	$Qry = $List[$i];                                                        

	$InputFile = $ORFeomesPath ."/". $Qry . $ext;                             
	$Db = $BlastPath ."/". $Qry;                                          

	$cmd = `makeblastdb -in $InputFile -dbtype nucl -parse_seqids -out $Db`;  

	Progress($n, $i);
}
exit;