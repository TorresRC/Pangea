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
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

my ($Usage, $ProjectName, $List, $TrustedORFeome, $eVal, $PIdent, $CPUs, $Help, $PanGenome,
    $CoreGenome, $Bolean);

GetOptions(
           'help'        => \$Help,
           'project|p=s' => \$ProjectName,
           'list|l=s'    => \$List,
           'trusted|c=s' => \$TrustedORFeome,
           'evalue|e=f'  => \$eVal,
           'ident|i=i'   => \$PIdent,
           'cpus|t=i'    => \$CPUs,
           'con-pan|P'   => \$PanGenome,
           'alncore|C'   => \$CoreGenome,
           'boleant|b'   => \$Bolean,
           ) or die "USAGE:\n  $0 [--help] [--project -p prefix] [--list -l filename]
      [--trusted -c filename] [--evalue -e evalue] [--ident -i integer]
      [--cpus -t integer]
\n  Use \'--help\' to print detailed descriptions of options.\n\n";

if($Help){
print "
\t--project <Project_Name>
\t--list <List_File_Name>
\t--trusted <TrustedFile.fasta>
\t--e-value <e-value>
\t--ident <Percetage_of_Identity>
\t--cpus <CPUs>\n\n";
exit;
}

my ($MainPath, $Project, $Src, $SortGenes, $FilterORFeomes, $MakeBlastDb, $InitialComparison, $GeneContent, $GeneContentPlot,
    $BoleanPresenceAbsence, $ConsensusPanGenome, $CoreAlign, $Start, $End, $RunTime);

$Start = time();
$MainPath = '/Users/rc/CoreGenome';
$Src      = $MainPath ."/". "src";

#$Script1  = $Src ."/". "1.-FormatFeatures.pl";
#$Script2  = $Src ."/". "2.-FilterFeatures.pl";
#$Script3  = $Src ."/". "3.-GetORFeomes.pl";

$SortGenes  = $Src ."/". "1.-SortGenes.pl";
$FilterORFeomes  = $Src ."/". "2.-FilterORFeomes.pl";
$MakeBlastDb  = $Src ."/". "4.-MakeBlastDBs.pl";
$InitialComparison  = $Src ."/". "5.-PreCoreGenome.pl";
$GeneContent  = $Src ."/". "6.-CoreGenome.pl";
$GeneContentPlot  = $Src ."/". "Plot.pl";
$BoleanPresenceAbsence = $Src ."/". "FormatPresenceAbsenceReport.pl";
$ConsensusPanGenome = $Src ."/". "ConsensusPanGenome.pl";
$CoreAlign = $Src ."/". "CoreAlign.pl";

#system("perl $SortGenes $ProjectName $List $TrustedORFeome $MainPath");
#system("perl $FilterORFeomes $ProjectName $List $TrustedORFeome $MainPath");
##system("perl $Script3 $ProjectName $List");
#system("perl $MakeBlastDb $ProjectName $List $TrustedORFeome $MainPath");
#system("perl $InitialComparison $ProjectName $List $TrustedORFeome $eVal $PIdent $CPUs $MainPath");
#system("perl $GeneContent $ProjectName $List $CPUs $MainPath");
system("perl $GeneContentPlot $ProjectName $List $MainPath");

if($PanGenome){
        system("perl $ConsensusPanGenome $ProjectName $List $MainPath");
}elsif($CoreGenome){
        system("perl $CoreAlign $ProjectName $List $MainPath");
}elsif($Bolean){
        system("perl $BoleanPresenceAbsence $ProjectName $List $MainPath");
}



$End = time();
$RunTime = ((($End - $Start)/60)/60);

print "\n\tFinished! The $ProjectName gene content analysis took $RunTime hours\n\n";
exit;