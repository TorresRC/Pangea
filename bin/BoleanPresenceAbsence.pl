#!/usr/bin/perl -w
#################################################################################
#Scipt FormatPresenceAbsence.pl                                                 #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          11 de octubre de 2017                                           #
#################################################################################
use strict; 
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $List, $MainPath);

$Usage = "\tUsage: CoreGenome.pl <Main_Path> <Project_Name>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$MainPath    = $ARGV[0];
$ProjectName = $ARGV[1];

my ($Project, $PresenceAbsenceFile, $BoleanReport, $TotalQry,
    $LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, $Row, $Field,
    $BoleanInformativeReport, $Line, $nElements, $Count, $ORF, $HeatMapRScript,
    $HeatMap, $Matrix);
my ($i, $j);
my (@PresenceAbsence, @PresenceAbsenceFields, @PresenceAbsenceArray,
    @Elements, @PresenceAbsenceMatrix);
my $BoleanTable = [ ];
my $GenesAnnotationReport = [ ];

$Project                 = $MainPath ."/". $ProjectName;
$PresenceAbsenceFile     = $Project ."/". $ProjectName . "_Presence_Absence.csv";
$BoleanReport            = $Project ."/". $ProjectName . "_Bolean_PanGenome.csv";
$BoleanInformativeReport = $Project ."/". $ProjectName . "_Bolean_AccessoryGenes.csv";
$HeatMapRScript          = $Project ."/". $ProjectName . "_Presence_Absence_HeatMapScript" . ".R";
$HeatMap                 = $Project ."/". $ProjectName . "_PanGenome_HeatMap" . ".pdf";

print "\nLoading the Presence/Absence file...";
($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, @PresenceAbsenceMatrix) = Matrix($PresenceAbsenceFile);
print "Done!";

print "\nProcessing";
$BoleanTable -> [0][0] = $PresenceAbsenceMatrix[0][0];
for ($i=0; $i<$LinesOnPresenceAbsence; $i++){
     $BoleanTable -> [$i][0] = $PresenceAbsenceMatrix[$i][0];
     for ($j=0; $j<$ColumnsOnPresenceAbsence; $j++){
          $BoleanTable -> [0][$j] = $PresenceAbsenceMatrix[0][$j];  
     }
}

#Setting the bolean data ("1" for presence and "0" for absence)
for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
     for ($j=1; $j<$ColumnsOnPresenceAbsence; $j++){
          if (defined $PresenceAbsenceMatrix[$i]->[$j]){          
               $Field = $PresenceAbsenceMatrix[$i]->[$j];
               if ($Field ne ""){   
                    $BoleanTable -> [$i][$j] = "1";
               }else{
                    $BoleanTable -> [$i][$j] = "0";
               }
          }else{
               $BoleanTable -> [$i][$j] = "0";
          }
     }
     Progress($LinesOnPresenceAbsence,$i);
}

print "Writing file...";
open (FILE, ">$BoleanReport");
for ($i=0; $i<$LinesOnPresenceAbsence; $i++){
     for ($j=0; $j<$ColumnsOnPresenceAbsence; $j++){
        if ($j < $ColumnsOnPresenceAbsence-1){
                print FILE $BoleanTable -> [$i][$j], ",";
        }else{
                print FILE $BoleanTable -> [$i][$j];
        }
     }
     print FILE "\n";
     Progress($LinesOnPresenceAbsence,$i);
}
close FILE;
print "Done!\n";

print "\nWriting a bolean file with only the informative ORFs...";
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

                $Count = 0;
                foreach my $Element(@Elements){
                        if ($Element ne "0"){
                                $Count++;
                        }
                }
                if($Count < $nElements && $Count > 1){
                        print FILE "$Line\n";   
                }
        }
close FILE;
print "Done!\n";


#Heatmap
open(RSCRIPT, ">$HeatMapRScript");
        print RSCRIPT 'library(gplots)' . "\n";
        print RSCRIPT 'library(RColorBrewer)' . "\n";
        
        print RSCRIPT "pdf(\"$HeatMap\", height=12, width=10)" . "\n";
        
        print RSCRIPT 'Colors <- colorRampPalette(c("beige","red"))(n=100)' . "\n";
        
        $Matrix = "Matrix";
        print RSCRIPT "df <- read.csv(\"$BoleanReport\")" . "\n";
        print RSCRIPT 'rnames <- df[,1]' . "\n";
        print RSCRIPT "$Matrix <- data.matrix(df[,2:ncol(df)])" . "\n";
        print RSCRIPT "rownames($Matrix) <- rnames" . "\n";
        print RSCRIPT "heatmap.2($Matrix," . "\n";
        print RSCRIPT "main = \"Pan-Genome of $ProjectName\"," . "\n";   
        print RSCRIPT 'density.info="none",' . "\n";                     
        print RSCRIPT 'trace = "none",' . "\n";                          
        print RSCRIPT 'xlab = "Class",' . "\n";
        print RSCRIPT 'ylab = "Genes",' . "\n";    
        print RSCRIPT 'Rowv = "NA",' . "\n";                             
        print RSCRIPT "Colv = as.dendrogram(hclust(as.dist(1-cor(t($Matrix)))))," . "\n";
        print RSCRIPT 'dendrogram = "column",' ."\n";                    
        
        print RSCRIPT 'adjCol= c(1,1),' ."\n";
        print RSCRIPT 'cexCol=0.3,' ."\n";                               
        print RSCRIPT 'cexRow=0.3,' ."\n";                               
        print RSCRIPT 'lhei=c(0.2,1),' ."\n";
        print RSCRIPT 'offsetRow=0,' ."\n";
        print RSCRIPT 'offsetCol=0,' ."\n";
        print RSCRIPT 'col = Colors)' . "\n";                            
        
        print RSCRIPT 'dev.off()';
close RSCRIPT;
system ("R CMD BATCH $HeatMapRScript");
system ("rm $HeatMapRScript");
        
system ("rm $Project/*.Rout $Project/Rplots.pdf");
print "Done!\n\n";

exit;