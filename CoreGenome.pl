#!/usr/bin/perl -w
#################################################################################
#   Programa CoreGenome                                                         #
#   Nota: Se debe ajustar la ruta de lib y de la variable $PathSeq a las que    #
#   realmente se tengan donde se instalción el programa.                        #
#                                                                               #
# Programador:   Alfonso Méndez Tenorio, Roberto C. Torres                      #
# Fecha:         2 de julio de 2017                                             #
#################################################################################
use strict;

my ($Usage, $ProjectName, $List, $TrustedORFeome, $eVal, $PIdent, $CPUs);

$Usage = "\tUsage: CoreGenome.pl <Project Name> <List File Name> <TrustedFile.fasta> <e-value> <Percentage of Identity> <CPUs>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName    = $ARGV[0];
$List           = $ARGV[1];
$TrustedORFeome = $ARGV[2];
$eVal           = $ARGV[3];
$PIdent         = $ARGV[4];
$CPUs           = $ARGV[5];

my ($MainPath, $Project, $Src, $Script1, $Script2, $Script3, $Script4, $Script5,
    $Script6, $Script7, $Start, $End, $RunTime);

$Start = time();
$MainPath = '/Users/rc/CoreGenome';
$Src      = $MainPath ."/". "src";
#$Script1  = $Src ."/". "1.-FormatFeatures.pl";
#$Script2  = $Src ."/". "2.-FilterFeatures.pl";
#$Script3  = $Src ."/". "3.-GetORFeomes.pl";
$Script1  = $Src ."/". "1.-SortGenes.pl";
$Script2  = $Src ."/". "2.-FilterORFeomes.pl";
$Script4  = $Src ."/". "4.-MakeBlastDBs.pl";
$Script5  = $Src ."/". "5.-PreCoreGenome.pl";
$Script6  = $Src ."/". "6.-CoreGenome.pl";
$Script7  = $Src ."/". "Plot.pl";

system("perl $Script1 $ProjectName $List $TrustedORFeome");
system("perl $Script2 $ProjectName $List $TrustedORFeome");
#system("perl $Script3 $ProjectName $List");
system("perl $Script4 $ProjectName $List $TrustedORFeome");
system("perl $Script5 $ProjectName $List $TrustedORFeome $eVal $PIdent $CPUs");
system("perl $Script6 $ProjectName $List $CPUs");
system("perl $Script7 $ProjectName $List");

$End = time();
$RunTime = ((($End - $Start)/60)/60);

print "\n\tFinished!, The $ProjectName gene content analysis took $RunTime hours\n\n";
exit;