#!/usr/bin/perl -w

#################################################################################
#By:       Roberto C. Torres & Mauricio Flores                                  #
#e-mail:   torres.roberto.c@gmail.com                                           #
#################################################################################
use strict;
use List::Util qw(max min);
use List::MoreUtils qw(any uniq first_index);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;
my $MainPath = "$FindBin::Bin";

my ($Usage, $TrainingFile, $MetadataFile, $OutPath, $Method, $X2, $IG, $OR,
    $PsCounts, $MI, $DotPlot, $HeatMapPlot, $Correlation, $Sort, $Clusters,
    $Dendrogram, $Runs, $Threshold, $AnnotatedPresenceAbsence);

$Usage = "\nUSAGE\n  $FindBin::Script <Observed Data [Absolute Path]>
                            <Metadata [Absolute Path]>
                            <Output Path [Relative Path]>
                            <Pseudo Counts Increase [Integer]>
                            <ChiSquare Test [Bolean]>
                            <Maximum Likelihood Estimation [Bolean]>\n\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$TrainingFile   = $ARGV[0];
$MetadataFile   = $ARGV[1];
$OutPath        = $ARGV[2];
$PsCounts       = $ARGV[3];
$Method         = $ARGV[4];  # <- X2, IG, OR, LO, MI
$DotPlot        = $ARGV[5];  # <- AllClasses, ForClass
$HeatMapPlot    = $ARGV[6];  # <- on/off
$Dendrogram     = $ARGV[7]; # <- off/row/column/both
$Threshold      = $ARGV[8];
$AnnotatedPresenceAbsence = $ARGV[9];

my($Test, $Run, $TestReport, $PercentagesReport, $Plot, $HeatMap, $PlotRScript,
   $LinesOnTrainingFile, $nFeature, $Line, $ColumnsOnTrainingFile, $N,
   $LinesOnMetaDataFile, $ColumnsOnMetaDataFile, $PossibleClass, $Column, $Class,
   $nClasses, $Element, $GlobalHits, $Hit, $Feature, $iClass, $a, $b, $c, $d,
   $nConfusion, $ChiConfidence, $Round, $HeatMapRScript, $Matrix, $CombinedInformative,
   $TestInformative, $PresenceInformative, $InformativeFeatures, $InformativeLines,
   $CombinedReport, $TestReportLine, $CombinedReportLine, $BoleanInformative,
   $TrainingHeader, $PresenceReportLine, $SummaryReport, $InformativeIndex,
   $InformativeClass, $SuperHeatMapRScript, $Percentage, $LinesOnAnnotation,
   $SelectedClass, $PresentSummaryReport, $AbsentSummaryReport, $BoleanInformativePresent,
   $BoleanInformativeAbsent, $TestInformativePresent, $TestInformativeAbsent);
my($i, $j, $nPresent);
my(@TrainingFile, @TrainingFileFields, @TrainingMatrix, @MetaDataFile,
   @MetaDataFileFields, @MetaDataMatrix, @Classes, @Elements, @ChiConfidence,
   @ChiConfidences, @TestReport, @Combined, @TestReportData, @Presence,
   @PresenceData, @Annotation, @AnnotationData, @Present, @Absent);
my(%ClassOfElement, %Elements, %pClass, %cpClass, %ClassHits,
   %HitsOfFeaturesInClass, %TotalFeatureHits, %Test,%PercentageOfFeatureInClass,
   %Gene, %EC, %Function, %InformativeClass, %Color, %InformativePresentClassCount,
   %InformativeAbsentClassCount, %ClassCount);
my(%a, %b, %c, %d);

my $Report = [ ];
my $Percentages = [ ];
my $Combined = [ ];

if ($Method eq "IG"){
   $Test = "Information_Gain";
}elsif ($Method eq "X2"){
   $Test = "Chi_Square";
}elsif ($Method eq "OR"){
   $Test = "Odds_Ratio";
}elsif ($Method eq "MI"){
   $Test = "Mutual_Information";
}elsif ($Method eq "LO"){
   $Test = "LogOdds";
}else {
   print "\nYou should select only one test option (--X2, --MLE or --OR)\n\tProgram finished!\n\n";
   exit;
}

# Loading the metadata file
($LinesOnMetaDataFile, $ColumnsOnMetaDataFile, @MetaDataMatrix) = Matrix($MetadataFile);

# Obtaining classes
print "\nThe following columns were detected as possible classes:";
for ($i=1;$i<$ColumnsOnMetaDataFile;$i++){
        $PossibleClass = $MetaDataMatrix[0][$i];
        print "\n\t[$i] $PossibleClass";
}
print "\n\nPlease type the index of the desired class: ";
$Column = <STDIN>;
chomp $Column;

$SelectedClass = $MetaDataMatrix[0]->[$Column];

for ($i=1;$i<$LinesOnMetaDataFile;$i++){
	$Class = $MetaDataMatrix[$i]->[$Column];
        #$Class =~ s/\s//g;
	push @Classes, $Class;
}

@Classes = uniq(@Classes);
$nClasses = scalar@Classes;
        
$TestReport               = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "Values" . ".csv";
$PercentagesReport        = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "PresencePercentages" . ".csv";
$CombinedReport           = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "ValuesAndPresencePercentages" . ".csv";
        
$TestInformative          = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "ValuesOfInformativeFeatures" . ".csv";
$TestInformativePresent   = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "ValuesOfInformativePresentFeatures" . ".csv";
$TestInformativeAbsent    = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "ValuesOfInformativeAbsentFeatures" . ".csv";
$PresenceInformative      = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "PercentagesOfInformativeFeatures" . ".csv";
$CombinedInformative      = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "ValuesAndPresencePercentagesOfInformativeFeatures" . ".csv";
$BoleanInformativePresent = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "BoleanInformativePresentFeatures" . ".csv";
$BoleanInformativeAbsent  = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "BoleanInformativeAbsentFeatures" . ".csv";
        
$SummaryReport            = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "AssociatedGenes" . ".csv";
$PresentSummaryReport     = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "AssociatedPresentGenes" . ".csv";
$AbsentSummaryReport      = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "AssociatedAbsentGenes" . ".csv";
        
$Plot                     = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "DotPlot" . ".pdf";
$HeatMap                  = $OutPath ."/". $SelectedClass ."_". $Test ."_". $Threshold . "CI" ."_". "HeatMap" . ".pdf";
$PlotRScript              = $OutPath ."/". $SelectedClass ."_". "DotPlotScript" . ".R";
$HeatMapRScript           = $OutPath ."/". $SelectedClass ."_". "HeatMapScript" . ".R";
$SuperHeatMapRScript      = $OutPath ."/". $SelectedClass ."_". "SuperHeatMapScript" . ".R";
        
# Loading the bolean training file
($LinesOnTrainingFile, $ColumnsOnTrainingFile, @TrainingMatrix) = Matrix($TrainingFile);
$nFeature = $LinesOnTrainingFile-1;
$N = $ColumnsOnTrainingFile-1;
        
$TrainingHeader = join ',', @{$TrainingMatrix[0]} [0..$N];
#chop $TrainingHeader;
#chop $TrainingHeader;
        
open (BOLEANPRESENT, ">$BoleanInformativePresent");
open (BOLEANABSENT, ">$BoleanInformativeAbsent");
   print BOLEANPRESENT "$TrainingHeader\n";
   print BOLEANABSENT "$TrainingHeader\n";
close BOLEANPRESENT;
close BOLEANABSENT;
        
for ($i=0;$i<$nClasses;$i++){
   for ($j=1;$j<$ColumnsOnTrainingFile;$j++){
      $Element = $TrainingMatrix[0]->[$j];
      $Class = $MetaDataMatrix[$j]->[$Column];
      $ClassOfElement{$Element} = $Class;
                        
      if($Class eq $Classes[$i]){
         $Elements{$Classes[$i]}++; #   <-------- Number of elements in each class
      }
   }
}

# Hits into the training matrix
$GlobalHits = 0;
for ($i=1; $i<$LinesOnTrainingFile; $i++){
   for ($j=1; $j<$ColumnsOnTrainingFile; $j++){
      $Hit = $TrainingMatrix[$i][$j];
      if ($Hit != 0){
         $GlobalHits++; #   <----------------------- Total of hits
      }
   }
}

#Hits of each class
foreach $Class(@Classes){
   for ($i=1;$i<$ColumnsOnTrainingFile; $i++){
      $Element = "$TrainingMatrix[0][$i]";
      if ($ClassOfElement{$Element} eq $Class){
         for ($j=1;$j<$LinesOnTrainingFile;$j++){
            $ClassHits{$Class} += $TrainingMatrix[$j][$i]+$PsCounts;  # <- Total of hits in class
         }
      }
   }
}

#Hits of each feature in each class
foreach $Class(@Classes){
   for ($i=1;$i<$LinesOnTrainingFile;$i++){
      $Feature = $TrainingMatrix[$i][0];
      $TotalFeatureHits{$Feature} = 0;
      for ($j=1;$j<$ColumnsOnTrainingFile; $j++){         
         $Element = $TrainingMatrix[0][$j];
         $TotalFeatureHits{$Feature} += $TrainingMatrix[$i][$j]+$PsCounts; # <- Total Feature Hits
         if ($ClassOfElement{$Element} eq $Class){
            $HitsOfFeaturesInClass{$Feature}{$Class} += $TrainingMatrix[$i][$j]+$PsCounts; # <- Total Feature Hits in Class
         }
      }
         $PercentageOfFeatureInClass{$Feature}{$Class} = ($HitsOfFeaturesInClass{$Feature}{$Class})/($Elements{$Class})*100; # <- Percentage of Feature Hits in Class
   }
}

# Statistic tests
$Report -> [0][0] = "Feature";
$Percentages -> [0][0] = "Feature";
$iClass = 1;
$InformativeLines = 0; 
for ($i=0; $i<$nClasses; $i++){
   $Class = $Classes[$i];
   $Report -> [0][$i+1] = $Class;
   $Percentages -> [0][$i+1] = $Class;
   for ($j=1;$j<$LinesOnTrainingFile;$j++){
      $Feature = $TrainingMatrix[$j][0];
        
      $a = (($HitsOfFeaturesInClass{$Feature}{$Class}))+0.001; # hits de sonda a en clase a
      $b = (($TotalFeatureHits{$Feature}-$HitsOfFeaturesInClass{$Feature}{$Class}))+0.001; # Hits de sonda a que no están en clase A
      $c = (($Elements{$Class}-$HitsOfFeaturesInClass{$Feature}{$Class}))+0.001; # Numero de mismaches en clase A (numero de ceros en clase A)
      $d = ((($N-$Elements{$Class})-($TotalFeatureHits{$Feature}-$HitsOfFeaturesInClass{$Feature}{$Class})))+0.001; # Numero de ceros fuera de A
      $nConfusion = $a+$b+$c+$d;
                    
      if ($Method eq "IG"){         # <------------------ Maximum Likelihood Estimation -- Information Gain
         $Test{$Feature} = ((-1*(($a+$c)/$nConfusion))*log10(($a+$c)/$nConfusion))+
         (($a/$nConfusion)*log10($a/($a+$b)))+
         (($c/$nConfusion)*log10($c/($c+$d)));
      }elsif ($Method eq "X2"){     # <------------------------------------ Chi square
         $Test{$Feature} = (($nConfusion*(($a*$d)-($b*$c))**2))/(($a+$c)*($a+$b)*($b+$d)*($c+$d));
      }elsif ($Method eq "OR"){     # <------------------------------------ Odds Ratios
         $Test{$Feature} = ($a*$d)/($b*$c);
      }elsif ($Method eq "LO"){     # <------------------------------------ LogOdds
         $Test{$Feature} = log2(($a*$d)/($b*$c));
      }elsif ($Method eq "MI"){     # <------------------------------------ Mutual information
         $Test{$Feature} = log10(($a*$nConfusion)/(($a+$b)*($a+$c)));
      }
              
      $Report -> [$j][0] = $Feature;
      $Report -> [$j][$iClass] = $Test{$Feature};
              
      $Percentages -> [$j][0] = $Feature;
      $Percentages -> [$j][$iClass] = $PercentageOfFeatureInClass{$Feature}{$Class};
   }
   $iClass++;
}
        
# Building output file
open (FILE, ">$TestReport");
open (PERCENTAGES, ">$PercentagesReport");
   for ($i=0;$i<$LinesOnTrainingFile;$i++){
      for ($j=0;$j<$nClasses+1;$j++){
         if($j < $nClasses){
            print FILE $Report -> [$i][$j], ",";
            print PERCENTAGES $Percentages -> [$i][$j], ",";
         }elsif($j == $nClasses){
            print FILE $Report -> [$i][$j];
            print PERCENTAGES $Percentages -> [$i][$j];
         }
      }
      print FILE "\n";
      print PERCENTAGES "\n";
   }
close FILE;
close PERCENTAGES;
        
for ($i=0;$i<$LinesOnTrainingFile;$i++){
   for ($j=0;$j<$nClasses*2+1;$j++){
      if ($j<$nClasses+1){
         $Combined -> [$i][$j] = $Report -> [$i][$j];
      }else{
         $Combined -> [$i][$j] = $Percentages -> [$i][$j-$nClasses];
      }
   }
}
        
open (FILE, ">$CombinedReport");
for ($i=0;$i<$LinesOnTrainingFile;$i++){
   for ($j=0;$j<$nClasses*2+1;$j++){
      print FILE $Combined -> [$i][$j], ",";
   }
   print FILE "\n";
}
close FILE;
        
@TestReport = ReadFile($TestReport);
@Presence   = ReadFile($PercentagesReport);
@Combined   = ReadFile($CombinedReport);
@Annotation = ReadFile($AnnotatedPresenceAbsence);
$LinesOnAnnotation = scalar@Annotation;
        
for ($i=1;$i<$LinesOnAnnotation;$i++){
   @AnnotationData = split(",",$Annotation[$i]);
   chomp@AnnotationData;
   $Feature            = $AnnotationData[0];
   chomp$Feature;
   $Gene{$Feature}     = $AnnotationData[1];
   $EC{$Feature}       = $AnnotationData[2];
   $Function{$Feature} = $AnnotationData[3];
}
        
for ($i=0;$i<$LinesOnTrainingFile;$i++){
   $TestReportLine     = $TestReport[$i];
   $PresenceReportLine = $Presence[$i];
   $CombinedReportLine = $Combined[$i];
                
   @PresenceData = split(",",$PresenceReportLine);
   shift@PresenceData;
                
   if ($i == 0){
      open (CHI, ">$TestInformative");
      open (CHIPRESENT, ">$TestInformativePresent");
      open (CHIABSENT, ">$TestInformativeAbsent");
      open (PRESENCE, ">$PresenceInformative");
      open (COMBINED, ">$CombinedInformative");
      open (PRESENT, ">$PresentSummaryReport");
      open (ABSENT, ">$AbsentSummaryReport");
         print CHIPRESENT $TestReportLine, "\n";
         print CHIABSENT $TestReportLine, "\n";
         print CHI $TestReportLine, "\n";
         print PRESENCE $PresenceReportLine, "\n";
         print COMBINED $CombinedReportLine, "\n";
         print PRESENT "ORF,Gene,EC_Number,Product,";
         print ABSENT "ORF,Gene,EC_Number,Product,";
         foreach $Class(@Classes){
            #chop$Class;
            print PRESENT "$Class(%),";
            print ABSENT "$Class(%),";
         }
         print PRESENT "Represented_Class,Statistical_Prediction\n";
         print ABSENT "Represented_Class,Statistical_Prediction\n";
      close CHIPRESENT;
      close CHIABSENT;
      close PRESENCE;
      close COMBINED;
      close PRESENT;
      close ABSENT;
      close CHI;
   }else{
      @TestReportData = split(",",$TestReportLine);
      $Feature = $TestReportData[0];
      chomp$Feature;
      shift@TestReportData;
                        
      if ($Method eq "X2"){
         if ( any { $_ > $Threshold} @TestReportData){
            open (CHI, ">>$TestInformative");                     
            open (PRESENCE, ">>$PresenceInformative");
            open (COMBINED, ">>$CombinedInformative");
               print CHI $TestReportLine, "\n";
               print PRESENCE $PresenceReportLine, "\n";
               print COMBINED $CombinedReportLine, "\n";
            close PRESENCE;
            close COMBINED;
            close CHI;
                                        
            $InformativeIndex = first_index{$_ eq max@TestReportData} @TestReportData;
            $InformativeClass = $Classes[$InformativeIndex];
            $InformativeClass{$Feature} = $InformativeClass;
            #open (SUM, ">>$SummaryReport");
            open (PRESENT, ">>$PresentSummaryReport");
            open (ABSENT, ">>$AbsentSummaryReport");
            $Function{$Feature} =~ s/\s//g;
            
            if ($PresenceData[$InformativeIndex] == max@PresenceData){
               push @Present, "$Feature";
               $InformativePresentClassCount{$Classes[$InformativeIndex]} += 1; 

               print PRESENT "$Feature,$Gene{$Feature},$EC{$Feature},$Function{$Feature},";
               foreach $Percentage(@PresenceData){
                  print PRESENT "$Percentage,";
               }
               print PRESENT "$Classes[$InformativeIndex],Present\n";
               
               open (CHIPRESENT, ">>$TestInformativePresent");
                  print CHIPRESENT $TestReportLine, "\n";
               close CHIPRESENT;
               
               open (BOLEAN, ">>$BoleanInformativePresent");
                  my $InformativeFeature = `grep -w $Feature $TrainingFile`;
                  chop $InformativeFeature;
                  chop $InformativeFeature;
                  print BOLEAN $InformativeFeature, "\n";
               close BOLEAN;
               
            }else{
               push @Absent, "$Feature";
               $InformativeAbsentClassCount{$Classes[$InformativeIndex]} += 1; 
               #print SUM "$Feature,$Classes[$InformativeIndex],$PresenceData[$InformativeIndex],Absent\n";                                 
               print ABSENT "$Feature,$Gene{$Feature},$EC{$Feature},$Function{$Feature},";
               foreach $Percentage(@PresenceData){
                  print ABSENT "$Percentage,";
               }
               print ABSENT "$Classes[$InformativeIndex],Absent\n";
               
               open (CHIABSENT, ">>$TestInformativeAbsent");
                  print CHIABSENT $TestReportLine, "\n";
               close CHIABSENT;
               
               open (BOLEAN, ">>$BoleanInformativeAbsent");
                  my $InformativeFeature = `grep -w $Feature $TrainingFile`;
                  chop $InformativeFeature;
                  chop $InformativeFeature;
                  print BOLEAN $InformativeFeature, "\n";
               close BOLEAN;
               
            }
            close PRESENT;
            close ABSENT;
         }
      }elsif($Method eq "LO"){
         if ( any { $_ > "3" or $_ < "-3"} @TestReportData){
            open (LO, ">>$TestInformative");
            open (PRESENCE, ">>$PresenceInformative");
            open (COMBINED, ">>$CombinedInformative");
               print LO $TestReportLine, "\n";
               print PRESENCE $PresenceReportLine, "\n";
               print COMBINED $CombinedReportLine, "\n";
            close LO;
            close PRESENCE;
            close COMBINED;
                                           
            open (BOLEAN, ">>$BoleanInformative");
               my $InformativeFeature = `grep -w $Feature $TrainingFile`;
               chop $InformativeFeature;
               chop $InformativeFeature;
               print BOLEAN $InformativeFeature, "\n";
            close BOLEAN;
         }
      }       
   }    
}

$nPresent = scalar@Present;
        
        # Building dot plot
        print "\n Building Plots...";
        chdir($OutPath);
        
        # Dot plot all clases
        if ($DotPlot eq "AllClasses"){
                open(RSCRIPT, ">$PlotRScript");
                        print RSCRIPT 'library(ggplot2)' . "\n";
                        print RSCRIPT "df <- read.csv(\"$TestReport\")" . "\n";
                        print RSCRIPT 'ggplot(df, aes(Feature))';
                        foreach $Class(@Classes){
                                print RSCRIPT "+ geom_point(aes(y=$Class,color=\"$Class\"))";
                        }
                        if($Method eq "X2"){
                                @ChiConfidences = (0.9,0.95,0.995,0.999);
                                foreach $ChiConfidence(@ChiConfidences){
                                        print RSCRIPT "+ geom_hline(aes(yintercept = qchisq($ChiConfidence, df=$nClasses-1), linetype=\"$ChiConfidence\"))";
                                }
                        }
                        print RSCRIPT "+ labs(x=\"Features\", y=\"Chi Values\", title= \"$Test\", color=\"Class\", linetype=\"Confidence Intervals\")";
                        if($N > 100){
                                print RSCRIPT '+ theme(axis.text.x = element_text(angle = 90, size=4, hjust = 1))' . "\n";
                        }
                        print RSCRIPT "\n";
                        print RSCRIPT "ggsave(\"$Plot\")" . "\n";
                        
                close RSCRIPT;
                system ("R CMD BATCH $PlotRScript");
                system ("rm $PlotRScript");
        }
        
        # Dot plot for class
        if ($DotPlot eq "ForClass"){
                foreach $Class(@Classes){
                        my $ClassPlotRScript = $OutPath ."/". $Class . "_DotPlotScript.R";
                        my $ClassPlot = $OutPath ."/". $Test ."_". $Class . "_DotPlot.pdf";
                        open(FILE, ">$ClassPlotRScript");
                                print FILE 'library(ggplot2)' . "\n";
                                print FILE "df <- read.csv(\"$TestReport\")" . "\n";
                                print FILE 'ggplot(df, aes(Feature))';
                                print FILE "+ geom_point(aes(y=$Class,color=\"$Class\"))";
                                if($Method eq "X2"){
                                        @ChiConfidences = (0.9,0.95,0.975,0.99,0.999);
                                        foreach $ChiConfidence(@ChiConfidences){
                                                print FILE "+ geom_hline(aes(yintercept = qchisq($ChiConfidence, df=$nClasses-1), linetype=\"$ChiConfidence\"))";
                                        }
                                }
                                print FILE "+ labs(x=\"Features\", y=\"Chi_Values\", title= \"$Test\", color=\"Class\", linetype=\"Confidence Intervals\")";
                                if($N > 100){
                                        print FILE '+ theme(axis.text.x = element_text(angle = 90, size=4, hjust = 1))' . "\n";
                                }
                                print FILE "\n";
                                print FILE "ggsave(\"$ClassPlot\")" . "\n";
                                
                        close FILE;
                        system ("R CMD BATCH $ClassPlotRScript");
                        system ("rm $ClassPlotRScript");
                }
        }
         
        # Heat Map;
        if ($HeatMapPlot eq "on"){
                open(RSCRIPT, ">$HeatMapRScript");
                        print RSCRIPT 'library(gplots)' . "\n";
                        print RSCRIPT 'library(RColorBrewer)' . "\n"; 
                        
                        print RSCRIPT 'Features <- c(';
                        for( $i=0; $i<$nPresent-1; $i++){
                           $Feature = $Present[$i];
                              print RSCRIPT "\"$Feature\",";
                        }
                        print RSCRIPT "\"$Present[$#Present]\")" . "\n";
                        
                        foreach $Class(@Classes){
                           print RSCRIPT "$Class <- c(";
                           for ($i=0; $i<$nPresent; $i++){
                              $Feature = $Present[$i];
                              if($InformativeClass{$Feature} eq $Class){
                                 $ClassCount{$Class} += 1;
                                 
                                 if($ClassCount{$Class} < $InformativePresentClassCount{$Class}){
                                    print RSCRIPT "\"$Feature\",";
                                 }else{
                                    print RSCRIPT "\"$Feature\")" . "\n";
                                 }
                              }
                           }
                        }
                        
                        print RSCRIPT "ClassesCol <- rep(\"black\",$nPresent)" . "\n";
                        
                        ####################
                        $Color{None_Atrophic_Gastritis} = "darkgreen";
                        $Color{Atrophic_Gastritis} = "gold";
                        $Color{Intestinal_Metaplasia} = "orange";
                        $Color{Gastric_Cancer} = "red";
                        ####################
                        foreach $Class(@Classes){
                           print RSCRIPT "ClassesCol[Features %in% $Class] <- \"$Color{$Class}\"" . "\n";
                        }
                        
                        print RSCRIPT "pdf(\"$HeatMap\",height=9,width=7)" . "\n";
                        
                        #print RSCRIPT 'Colors <- colorRampPalette(c("red", "yellow", "green"))(n=100)' . "\n";
                        #print RSCRIPT 'Colors <- colorRampPalette(c("oldlace", "blue"))(n=150)' . "\n";
                        print RSCRIPT 'Colors <- colorRampPalette(c("oldlace", "red"))(n=150)' . "\n";
                        
                        $Matrix = "Matrix";
                        print RSCRIPT "df <- read.csv(\"$TestInformativePresent\")" . "\n";
                        print RSCRIPT 'rnames <- df[,1]' . "\n";
                        print RSCRIPT "$Matrix <- data.matrix(df[,2:ncol(df)])" . "\n";
                        print RSCRIPT "rownames($Matrix) <- rnames" . "\n";
                        print RSCRIPT "hc <- hclust(as.dist(1-cor(t($Matrix))))" . "\n";
                        
                        print RSCRIPT "heatmap.2($Matrix," . "\n";
                        print RSCRIPT "Rowv=as.dendrogram(hc)," . "\n";
                        print RSCRIPT "main = \"$Test\"," . "\n";                        # Title
                        print RSCRIPT 'key.title = NA,' . "\n";
                        print RSCRIPT "key.xlab = \"$Test Values\"," . "\n";
                        print RSCRIPT 'trace = "none",' . "\n";                          # Turns of trace lines in heat map
                        print RSCRIPT 'tracecol = "red",' . "\n";
                        print RSCRIPT 'srtCol= 0,' ."\n";
                        print RSCRIPT 'adjCol= c (0.5,1),' ."\n";
                        print RSCRIPT 'cexCol=0.8,' ."\n";
                        print RSCRIPT 'lhei=c(0.2,1),' ."\n";                        
                        
                        #$Round = 5;
                        #if ($Correlation eq "on"){
                        #        print RSCRIPT 'xlab = "Class",' . "\n";
                        #        print RSCRIPT 'ylab = "Class",' . "\n";
                        #        print RSCRIPT "cellnote = round($Matrix,$Round)," . "\n"; # Shows data in cell
                        #        print RSCRIPT 'srtRow= 90,' ."\n";
                        #        print RSCRIPT 'adjRow= c (0.5,1),' ."\n";
                        #        print RSCRIPT 'cexRow=0.5,' ."\n";
                        #}else{
                        #        print RSCRIPT 'xlab = "Class",' . "\n";
                        #        print RSCRIPT 'ylab = "Feature",' . "\n";
                        #        if ($nClasses < 11 && $nFeature < 101){
                        #                print RSCRIPT "cellnote = round($Matrix,$Round)," . "\n"; # Shows data in cell
                        #        }
                        #}
                        
                        if ($Dendrogram eq "off"){
                                print RSCRIPT 'dendrogram = "none",' ."\n";              # Hides dendrogram
                        }elsif($Dendrogram eq "row"){
                                print RSCRIPT 'dendrogram = "row",' ."\n";               # Shows dendrogram for rows
                        }elsif($Dendrogram eq "column"){
                                print RSCRIPT 'dendrogram = "column",' ."\n";            # Shows dendrogram for columns
                        }elsif($Dendrogram eq "both"){
                                print RSCRIPT 'dendrogram = "both",' ."\n";              # Shos dendrogram for rows and columns
                        }
                        
                        print RSCRIPT 'RowSideColor=ClassesCol,' ."\n";
                        #print RSCRIPT 'col=brewer.pal(4,"Blues"))' . "\n";               # Use defined palette
                        print RSCRIPT 'col=Colors)' . "\n";
                        
                        print RSCRIPT 'dev.off()';
                close RSCRIPT;
                system ("R CMD BATCH $HeatMapRScript");
                system ("rm $HeatMapRScript");
        }
           
        system ("rm $OutPath/*.Rout $OutPath/Rplots.pdf");
        
        print "Done!\n\n";
        
        $TrainingFile = $BoleanInformative;
        
exit;