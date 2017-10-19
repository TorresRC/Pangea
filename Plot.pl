#!/usr/bin/perl -w
#################################################################################
#   Programa CoreGenome                                                         #
#   Nota: Se debe ajustar la ruta de lib y de la variable $PathSeq a las que    #
#   realmente se tengan donde se instalción el programa.                        #
#                                                                               #
# Programador:   Roberto C. Torres                                              #
# Fecha:         18 de Octubre de 2017                                          #
#################################################################################
use strict;
use lib '/Users/rc/lib';
use Routines;

my ($Usage, $ProjectName, $List);

$Usage = "\tUsage: CoreGenome.pl <Project Name> <List File Name>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List = $ARGV[1];

my ($MainPath, $Project, $MainList, $PresenceAbsence, $Plot, $RScript);
my ($n);
my (@List);

$MainPath        = "/Users/rc/CoreGenome";
$Project         = $MainPath ."/". $ProjectName;
$MainList        = $Project ."/". $List;
$PresenceAbsence = $Project ."/". $ProjectName . '_Statistics.csv';
$Plot            = $Project ."/". $ProjectName . "_GeneContentPlot.pdf";
$RScript         = $Project ."/". "GeneContentScript.R";

@List = ReadFile($MainList);
$n = scalar@List;

chdir ($Project);

open (RSCRIPT, ">$RScript");
        print RSCRIPT 'library(ggplot2)' . "\n";
        print RSCRIPT "df <- read.csv(\"$PresenceAbsence\")" . "\n";
        print RSCRIPT 'ggplot(df, aes(NumberOfNewStrains))';
        print RSCRIPT '+ geom_line(aes(y=CoreGenome,linetype="CoreGenome"))';
        print RSCRIPT '+ geom_line(aes(y=PanGenome,linetype="PanGenome"))';
        print RSCRIPT '+ geom_line(aes(y=NewGenes,linetype="NewGenes"))';
        print RSCRIPT "+ scale_x_continuous(breaks = 0:$n+1)";
        print RSCRIPT "+ labs(x=\"Number of Genomes\", y=\"Number of Genes\", title= \"$ProjectName Gene Content\")";
        print RSCRIPT '+ scale_linetype_discrete(name=NULL)';
        print RSCRIPT '+ theme(axis.text.x = element_text(angle = 90, size=4, hjust = 1))' . "\n";
#        print RSCRIPT "ggsave(\"$Plot\")";
close RSCRIPT;

print "\n Building a Gene Content Plot...";
system ("R CMD BATCH $RScript");
print "Done!\n\n";

system ("rm $RScript $Project/*.Rout");

exit;