#!/usr/bin/perl -w

#################################################################################
#   Programa Filter Features                                                    #
#   Nota: Se debe ajustar la ruta de lib y de la variable $PathSeq a las que    #
#   realmente se tengan donde se instalción el programa.                        #
#                                                                               #
# Programador:   Alfonso Méndez Tenorio, Roberto C. Torres                      #
# Fecha:         10 de enero de 2017                                            #
#################################################################################

use strict; 
use lib '/home/bioinformatica/CoreGenome/src/lib';
#use lib '/Users/Roberto/Documents/lib';
use rutinas;


my ($Usage, $ProjectName, $List);

$Usage = "\tUsage: GetORFeomes.pl <Project Name> <List File Name>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List = $ARGV[1];

my ($MainPath, $Project, $GenomePath, $FeaturesPath, $GeneExtractor,
    $InPath, $OutPath, $MainList, $Qry, $Prefix, $Genome, $Features, $SeqExt,
    $FeatureExt, $ORFeomeExt, $Out, $cmd, $LogFile);
my (@List, @Split);

#$MainPath = "/Users/Roberto/CoreGenome";
$MainPath = "/home/bioinformatica/CoreGenome";
$Project = $MainPath ."/". $ProjectName;
#$GenomePath     = $Project ."/". "Genomes";
$GenomePath     = $MainPath ."/". "Genomes";
#$FeaturesPath   = $Project ."/". "Features/ptt";
$FeaturesPath   = $MainPath ."/". "Features/ptt";
#$GeneExtractor  = $MainPath ."/". "src/GeneExtractorMacOS/gene_extractor";
$GeneExtractor  = $MainPath ."/". "src/gene_extractor";

$InPath     = $FeaturesPath ."/". "Filtered";
#$OutPath    = $Project ."/". 'ORFeomes';
$OutPath    = $MainPath ."/". 'ORFeomes';
$MainList   = $Project ."/". $List;
$SeqExt     = ".fasta";
$FeatureExt = ".ptt";
$ORFeomeExt = ".fasta";
$LogFile = $Project ."/". $ProjectName . ".log";


open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

MakeDir($OutPath);

@List = ReadFile($MainList);

foreach $Qry (@List){
	$Genome		= $GenomePath ."/". $Qry . $SeqExt;
	$Features	= $InPath ."/". "Filtered_" . $Qry . $FeatureExt;
	$Out		= $OutPath ."/". $Qry . $ORFeomeExt;
    
    print "\nGetting ORFeome of $Qry...";
    
	$cmd =`$GeneExtractor $Genome $Features $Out`;
    
    print "Done.";
}
close STDERR;
print "\n\n\tFinished!\n\n";
exit;
