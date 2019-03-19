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
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProgressFile, $Prefix, $OutPath, $List);

$Usage = "\tUsage: CoreGenome.pl <Main_Path> <Project_Name> <List_File_Name>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProgressFile = $ARGV[0];
$Prefix       = $ARGV[1];
$List         = $ARGV[2];
$OutPath      = $ARGV[3];

my ($Plot, $RScript);
my ($n);
my (@List);

$Plot     = $OutPath ."/". $Prefix . "_GeneContentPlot.pdf";
$RScript  = $OutPath ."/". $Prefix . "_GeneContentScript.R";

@List = ReadFile($List);
$n = scalar@List;

chdir ($OutPath);

open (RSCRIPT, ">$RScript");
        print RSCRIPT 'library(ggplot2)' . "\n";
        print RSCRIPT "df <- read.csv(\"$ProgressFile\")" . "\n";
        print RSCRIPT 'ggplot(df, aes(NumberOfNewStrains))';
        print RSCRIPT '+ geom_line(aes(y=CoreGenome,linetype="CoreGenome"))';
        print RSCRIPT '+ geom_line(aes(y=PanGenome,linetype="PanGenome"))';
        print RSCRIPT '+ geom_line(aes(y=NewGenes,linetype="NewGenes"))';
        print RSCRIPT "+ scale_x_continuous(breaks = 0:$n+1)";
        print RSCRIPT "+ labs(x=\"Number of Genomes\", y=\"Number of Genes\", title= \"$Prefix Gene Content\")";
        print RSCRIPT '+ scale_linetype_discrete(name=NULL)';
        print RSCRIPT '+ theme(axis.text.x = element_text(angle = 90, size=4, hjust = 1))' . "\n";
        print RSCRIPT "ggsave(\"$Plot\")";
close RSCRIPT;

print "\n Building a Gene Content Plot...";
system ("R CMD BATCH $RScript");
print "Done!\n\n";

system ("rm $RScript $OutPath/*.Rout $OutPath/*Rplots.pdf");

exit;