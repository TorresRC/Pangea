#!/usr/bin/perl -w
#################################################################################
#Scipt FormatPresenceAbsence.pl                                                 #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          11 de octubre de 2017                                           #
#################################################################################
use strict;
use List::Util qw(sum sum0);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $PresenceAbsenceFile, $OutPath);

$Usage = "\tUsage: CoreGenome.pl <Main_Path> <Project_Name>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName         = $ARGV[0];
$PresenceAbsenceFile = $ARGV[1];
$OutPath             = $ARGV[2];

my ($Project, $BoleanReport, $TotalQry, $LinesOnPresenceAbsence,
    $ColumnsOnPresenceAbsence, $Row, $Field, $BoleanInformativeReport, $Line,
    $nElements, $Count, $ORF, $HeatMapRScript, $HeatMap, $Matrix, $Sum, $Strain);
my ($i, $j, $N);
my (@PresenceAbsence, @PresenceAbsenceFields, @PresenceAbsenceArray,
    @Elements, @PresenceAbsenceMatrix);
my (%Count);
my $BoleanTable = [ ];
my $GenesAnnotationReport = [ ];


$BoleanReport            = $OutPath ."/". $ProjectName . "_Boolean_PanGenome.csv";
$BoleanInformativeReport = $OutPath ."/". $ProjectName . "_Boolean_AccessoryGenes.csv";
$HeatMapRScript          = $OutPath ."/". $ProjectName . "_Presence_Absence_HeatMapScript" . ".R";
$HeatMap                 = $OutPath ."/". $ProjectName . "_PanGenome_HeatMap" . ".pdf";

print "\nLoading the Presence/Absence file...";
($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, @PresenceAbsenceMatrix) = Matrix($PresenceAbsenceFile);
print "Done!";

$N = $ColumnsOnPresenceAbsence-1;

print "\nTransforming into binary data...";
$BoleanTable -> [0][0] = $PresenceAbsenceMatrix[0][0];
for ($i=0; $i<$LinesOnPresenceAbsence; $i++){
     $BoleanTable -> [$i][0] = $PresenceAbsenceMatrix[$i][0];
     for ($j=0; $j<$ColumnsOnPresenceAbsence; $j++){
          $BoleanTable -> [0][$j] = $PresenceAbsenceMatrix[0][$j];  
     }
}

for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
    $ORF = $PresenceAbsence[$i][0];
     for ($j=1; $j<$ColumnsOnPresenceAbsence; $j++){
        $Strain = $PresenceAbsence[0][$j];
        $Count{$ORF}{$Strain} = 0;
     }
}

#Setting the bolean data ("1" for presence and "0" for absence)
for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
    $ORF = $PresenceAbsence[$i][0];
     for ($j=1; $j<$ColumnsOnPresenceAbsence; $j++){
        $Strain = $PresenceAbsence[0][$j];
          if (defined $PresenceAbsenceMatrix[$i]->[$j]){          
               $Field = $PresenceAbsenceMatrix[$i]->[$j];
               $Field =~s/\W//g;
               if ($Field ne ""){
               #if ($Field eq 1){ 
                    $BoleanTable -> [$i][$j] = "1";
                    $Count{$ORF}{$Strain} += 1; 
               }else{
               #}elsif($Field eq 0){
                    $BoleanTable -> [$i][$j] = "0";
               }
          }else{
               $BoleanTable -> [$i][$j] = "0";
          }
     }
}
print "Done!";

print "\nBuilding a Presence/Absence binary file for the entire Pan-Genome...";
open (FILE, ">$BoleanReport");
for ($i=0; $i<$LinesOnPresenceAbsence; $i++){
    $ORF = $PresenceAbsence[$i][0];
     for ($j=0; $j<$ColumnsOnPresenceAbsence; $j++){
        $Strain = $PresenceAbsence[0][$j];
        if ($Count{$ORF}{$Strain} > 0){
            if ($j < $ColumnsOnPresenceAbsence-1){
                print FILE $BoleanTable -> [$i][$j], ",";
            }else{
                print FILE $BoleanTable -> [$i][$j];
            }
        }
     }
     print FILE "\n";
}
close FILE;
print "Done!";

print "\nBuilding a Presence/Absence binary file for accessory genes...";
my @BoleanReport = ReadFile($BoleanReport);  
open (FILE,">$BoleanInformativeReport");
        print FILE "$BoleanReport[0]\n";
        for ($i=1; $i<scalar@BoleanReport;$i++){
                $Line = $BoleanReport[$i];
                @Elements = split(",",$Line);
                chomp@Elements;
                $ORF = $Elements[0];
                shift@Elements;
                $nElements = scalar@Elements;
                
                my $SpecificCutoff = $nElements*0.05;
                my $CoreCutoff = $nElements*0.9;
                #if(sum(@Elements) < $nElements && sum(@Elements) > 1){
                if(sum(@Elements) < $CoreCutoff && sum(@Elements) > $SpecificCutoff){
                #if(sum(@Elements) > $AccessoryCutoff){
                        print FILE "$Line\n";
                }
        }
close FILE;
print "Done!";

#Heatmap
print "\nGenerating a plot for the Present/Absent Pan-Genome genes...";
open(RSCRIPT, ">$HeatMapRScript");
        print RSCRIPT 'library(gplots)' . "\n";
        print RSCRIPT 'library(RColorBrewer)' . "\n";
        
        print RSCRIPT "pdf(\"$HeatMap\", height=12, width=10)" . "\n";
        
        print RSCRIPT 'Colors <- colorRampPalette(c("beige","red"))(n=2)' . "\n";
        
        $Matrix = "Matrix";
        print RSCRIPT "df <- read.csv(\"$BoleanReport\")" . "\n";
        print RSCRIPT 'rnames <- df[,1]' . "\n";
        print RSCRIPT "$Matrix <- data.matrix(df[,2:ncol(df)])" . "\n";
        print RSCRIPT "rownames($Matrix) <- rnames" . "\n";
        print RSCRIPT "heatmap.2($Matrix," . "\n";
        print RSCRIPT "main = \"Pan-Genome of $ProjectName\"," . "\n";   
        print RSCRIPT 'density.info="none",' . "\n";
        print RSCRIPT 'key = FALSE,' . "\n"; 
        print RSCRIPT 'trace = "none",' . "\n";                          
        print RSCRIPT 'xlab = "Genomes",' . "\n";
        print RSCRIPT 'ylab = "Genes",' . "\n";    
        print RSCRIPT 'Rowv = "NA",' . "\n";                             
        print RSCRIPT 'Colv = "NA",' . "\n";
        print RSCRIPT 'dendrogram = "none",' ."\n";                    
        
        print RSCRIPT "labCol= c(1:ncol($Matrix))," ."\n";
        print RSCRIPT "labRow= c(1:nrow($Matrix))," ."\n";
        print RSCRIPT 'adjCol= c(1,1),' ."\n";
        print RSCRIPT 'cexCol=0.3,' ."\n";                               
        print RSCRIPT 'cexRow=0.3,' ."\n";
        print RSCRIPT 'lwid=c(1,6),' ."\n";
        print RSCRIPT 'lhei=c(1.5,10),' ."\n";
        print RSCRIPT 'offsetRow=0,' ."\n";
        print RSCRIPT 'offsetCol=0,' ."\n";
        print RSCRIPT 'col = Colors)' . "\n";                            
        
        print RSCRIPT 'dev.off()';
close RSCRIPT;
system ("R CMD BATCH $HeatMapRScript");
system ("rm $HeatMapRScript");
        
#system ("rm *.Rout");
print "Done!\n\n";
exit;