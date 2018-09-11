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
    $Dendrogram, $Runs, $Threshold, $AnnotatedPresenceAbsence, $ProjectName,
    $MappingFile);

$Usage = "\nUSAGE\n  $FindBin::Script <Observed Data [Absolute Path]>
                            <Metadata [Absolute Path]>
                            <Output Path [Relative Path]>
                            <Pseudo Counts Increase [Integer]>
                            <ChiSquare Test [Boolean]>
                            <Maximum Likelihood Estimation [Boolean]>\n\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;

$ProjectName              = $ARGV[0];
$TrainingFile             = $ARGV[1];
$MetadataFile             = $ARGV[2];
$OutPath                  = $ARGV[3];
$AnnotatedPresenceAbsence = $ARGV[4];
$Threshold                = $ARGV[5];
$MappingFile              = $ARGV[6];

#$Method                   = $ARGV[7];  # <- X2, IG, OR, LO, MI
#$DotPlot                  = $ARGV[8];  # <- AllClasses, ByClass

$Method = "X2";
$DotPlot = "AllClasses";


my($Test, $Run, $TestReport, $PercentagesReport, $Plot, $HeatMap, $PlotRScript,
   $LinesOnTrainingFile, $nFeature, $Line, $ColumnsOnTrainingFile, $N,
   $LinesOnMetaDataFile, $ColumnsOnMetaDataFile, $PossibleClass, $Column, $Class,
   $nClasses, $Element, $GlobalHits, $Hit, $Feature, $iClass, $a, $b, $c, $d,
   $nContingency, $ChiConfidence, $Round, $HeatMapRScript, $Matrix, $CombinedInformative,
   $TestInformative, $PresenceInformative, $InformativeFeatures, $InformativeLines,
   $CombinedReport, $TestReportLine, $CombinedReportLine, $BooleanAccessory,
   $TrainingHeader, $PresenceReportLine, $SignificantReport, $InformativeIndex,
   $InformativeClass, $SuperHeatMapRScript, $Percentage, $LinesOnAnnotation,
   $SelectedClass, $PresentSummaryReport, $AbsentSummaryReport, $BooleanAccessoryPresent,
   $BooleanAccessoryAbsent, $TestInformativePresent, $TestInformativeAbsent,
   $StatCondition, $ChiValue, $SignificantGene, $BooleanSignificantGene, $SignificantGeneData,
   $Strain, $nStrains, $PresentHeatMap, $PresentHeatMapRScript, $AbsentHeatMap,
   $AbsentHeatMapRScript, $CorrelatedHeatMap, $CorrelatedHeatMapRScript, $BooleanFile,
   $HeatMapOut, $ConsensusSeq, $Class1, $Class2, $Class3, $LinesOnTestReport,
   $ColumnsOnTestReport, $TestReportHeader, $NdClass, $RdClass, $nClasses2, $nClasses3,
   $PermutationsFile);
my($i, $j, $M, $nPresent, $nAbsent, $BooleanPresent, $BooleanAbsent, $nSignificant,
   $BooleanMostlyPresent, $BooleanMostlyAbsent);
my ($FILE, $PRESENT, $ABSENT, $CORRELATED, $BOOLEAN, $SIGNIFICANTPRESENT, $SIGNIFICANTABSENT);
my(@TrainingFile, @TrainingFileFields, @TrainingMatrix, @MetaDataFile,
   @MetaDataFileFields, @MetaDataMatrix, @Classes, @Elements, @ChiConfidence,
   @ChiConfidences, @TestReport, @Combined, @TestReportData, @Presence,
   @PresenceData, @Annotation, @AnnotationData, @Present, @Absent, @InformativeFiles,
   @Percentages, @ChiValues, @SignificantGenes, @Strains, @SignificantAbsent, @SignificantPresent,
   @Columns, @TestReportMatrix, @Classes2, @Classes3);
my(%ClassOfElement, %Elements, %pClass, %cpClass, %ClassHits,
   %HitsOfFeaturesInClass, %TotalFeatureHits, %Test,%PercentageOfFeatureInClass,
   %Gene, %EC, %Function, %InformativeClass, %Color, %PresentClassCount,
   %AbsentClassCount, %ClassCount, %PercentageOfHit, %Phenotypes, %Classes, %SignificantPresentClassCount,
   %SignificantAbsentClassCount, %ChiValues, %BooleanData, %TestReportData, %Classes2, %Classes3, %Percentages);
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

MakeDir($OutPath);

# Loading the metadata file
($LinesOnMetaDataFile, $ColumnsOnMetaDataFile, @MetaDataMatrix) = Matrix($MetadataFile);

# Obtaining classes
print "\nThe following columns were detected as possible classes:";
for ($i=1;$i<$ColumnsOnMetaDataFile;$i++){
        $PossibleClass = $MetaDataMatrix[0][$i];
        print "\n\t[$i] $PossibleClass";
}
print "\n\nPlease type the index of the desired class (1 to 3 classes could be selected using a comma separation but only the first one will be used in the clustering analysis. ej. 1,2,3): ";
$Column = <STDIN>;
chomp $Column;

@Columns = split(",",$Column);
$Class1 = $Columns[0];
if (scalar@Columns > 1){
   $Class2 = $Columns[1];
}
if (scalar@Columns > 2){
   $Class3 = $Columns[2];
}

$SelectedClass = $MetaDataMatrix[0]->[$Class1];


$TestReport               = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "ValuesOfAllGenes" . ".csv";
$SignificantReport        = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "ValuesOfSignificantGenes" . ".csv";
$BooleanAccessory         = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "BooleanSignificantGenes" . ".csv";
$BooleanMostlyPresent     = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "BooleanMostlyPresent" . ".csv";
$BooleanMostlyAbsent      = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "BooleanMostlyAbsent" . ".csv";
$ConsensusSeq             = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "SignificantGenes" . ".fasta";
       
$Plot                     = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "ValuesOfAllGenes_DotPlot" . ".pdf";
$PlotRScript              = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "DotPlotScript" . ".R";
$PresentHeatMapRScript    = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "PresentHeatMapScript" . ".R";
$AbsentHeatMapRScript     = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "AbsentHeatMapScript" . ".R";
$CorrelatedHeatMapRScript = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "CorrelatedHeatMapScript" . ".R";
$HeatMapOut               = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI";

for ($i=1;$i<$LinesOnMetaDataFile;$i++){
   $Strain = $MetaDataMatrix[$i]->[0];
   push @Strains, $Strain;
   
	$Class = $MetaDataMatrix[$i]->[$Class1];
   $ClassOfElement{$Strain} = $Class;
   $Elements{$Class}++; #   <---------------------------------------------------- Number of elements in each class
   push @Classes, $Class;
   push (@{$Classes{$Class}}, $Strain);
   
   if (scalar@Columns > 1){
      $NdClass = $MetaDataMatrix[$i]->[$Class2];
      push @Classes2, $NdClass;
      push (@{$Classes2{$NdClass}}, $Strain);
      
      $RdClass = $MetaDataMatrix[$i][$Class3];
      push @Classes3, $RdClass;
      push (@{$Classes3{$RdClass}}, $Strain);
   }   
}

@Classes = uniq(@Classes);
$nClasses = scalar@Classes;
$nStrains = scalar@Strains;

if (scalar@Columns > 1){
   @Classes2 = uniq(@Classes2);
   @Classes3 = uniq(@Classes3);
   $nClasses2 = scalar@Classes2;
   $nClasses3 = scalar@Classes3;
}

@TrainingFile = ReadFile($TrainingFile);
        
# Loading the boolean training file
($LinesOnTrainingFile, $ColumnsOnTrainingFile, @TrainingMatrix) = Matrix($TrainingFile);
$nFeature = $LinesOnTrainingFile-1;
$N = $ColumnsOnTrainingFile-1;
        
$TrainingHeader = join ',', @{$TrainingMatrix[0]} [0..$N];

for ($i=1; $i<$LinesOnTrainingFile; $i++){
   $Feature = $TrainingMatrix[$i][0];
   $BooleanData{$Feature} = join ',', @{$TrainingMatrix[$i]} [1..$N];
}
        
open (BOOLEANACCESSORY, ">$BooleanAccessory");
   print BOOLEANACCESSORY "$TrainingHeader\n";
close BOOLEANACCESSORY;
        
#for ($i=0;$i<$nClasses;$i++){
#   for ($j=1;$j<$ColumnsOnTrainingFile;$j++){
#      $Element = $TrainingMatrix[0]->[$j];
#      $Class = $MetaDataMatrix[$j]->[$Class1];
#      $ClassOfElement{$Element} = $Class;
#                        
#      if($Class eq $Classes[$i]){
#         $Elements{$Classes[$i]}++; #   <--------------------------------------- Number of elements in each class
#      }
#   }
#}

# Hits into the training matrix
$GlobalHits = 0;
for ($i=1; $i<$LinesOnTrainingFile; $i++){
   for ($j=1; $j<$ColumnsOnTrainingFile; $j++){
      $Hit = $TrainingMatrix[$i][$j];
      if ($Hit != 0){
         $GlobalHits++; #   <---------------------------------------------------- Total of hits
      }
   }
}

#Hits on each class
for ($i=1;$i<$ColumnsOnTrainingFile; $i++){
   $Element = "$TrainingMatrix[0][$i]";
   for ($j=1;$j<$LinesOnTrainingFile;$j++){
      $ClassHits{$ClassOfElement{$Element}} += $TrainingMatrix[$j][$i];  # <----- Total of hits in class
   }
}

#Hits of each feature in each class
for ($i=1;$i<$LinesOnTrainingFile;$i++){
   $Feature = $TrainingMatrix[$i][0];
   $TotalFeatureHits{$Feature} = 0;
   for ($j=1;$j<$ColumnsOnTrainingFile; $j++){
      $Element = $TrainingMatrix[0][$j];
      $TotalFeatureHits{$Feature} += $TrainingMatrix[$i][$j]; # <---------------- Total Feature Hits
      $HitsOfFeaturesInClass{$Feature}{$ClassOfElement{$Element}} += $TrainingMatrix[$i][$j];# <- Total Feature Hits in Class
   }
   $PercentageOfHit{$Feature} = ($TotalFeatureHits{$Feature}/$N)*100;
   foreach $Class(@Classes){
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
        
      $a = (($HitsOfFeaturesInClass{$Feature}{$Class}))+0.001; # <--------------- Hits de sonda a en clase a
      $b = (($TotalFeatureHits{$Feature}-$HitsOfFeaturesInClass{$Feature}{$Class}))+0.001; # <- Hits de sonda a que no están en clase A
      $c = (($Elements{$Class}-$HitsOfFeaturesInClass{$Feature}{$Class}))+0.001; # <----------- Numero de mismaches en clase A (numero de ceros en clase A)
      $d = ((($N-$Elements{$Class})-($TotalFeatureHits{$Feature}-$HitsOfFeaturesInClass{$Feature}{$Class})))+0.001; # <- Numero de ceros fuera de A
      #$nContingency = $a+$b+$c+$d;
      $nContingency = $N;
                    
      if ($Method eq "IG"){         # <------------------------------------------ Maximum Likelihood Estimation -- Information Gain
         $Test{$Feature} = ((-1*(($a+$c)/$nContingency))*log10(($a+$c)/$nContingency))+
         (($a/$nContingency)*log10($a/($a+$b)))+
         (($c/$nContingency)*log10($c/($c+$d)));
      }elsif ($Method eq "X2"){     # <------------------------------------------ Chi square
         $Test{$Feature} = (($nContingency*(($a*$d)-($b*$c))**2))/(($a+$c)*($a+$b)*($b+$d)*($c+$d));
      }elsif ($Method eq "OR"){     # <------------------------------------------ Odds Ratios
         $Test{$Feature} = ($a*$d)/($b*$c);
      }elsif ($Method eq "LO"){     # <------------------------------------------ LogOdds
         $Test{$Feature} = log2(($a*$d)/($b*$c));
      }elsif ($Method eq "MI"){     # <------------------------------------------ Mutual information
         $Test{$Feature} = log10(($a*$nContingency)/(($a+$b)*($a+$c)));
      }
              
      $Report -> [$j][0] = $Feature;
      $Report -> [$j][$iClass] = $Test{$Feature};
              
      $Percentages -> [$j][0] = $Feature;
      $Percentages -> [$j][$iClass] = $PercentageOfFeatureInClass{$Feature}{$Class};
   }
   $iClass++;
}

print "\nBuilding OutPut Files...";
@Annotation = ReadFile($AnnotatedPresenceAbsence);
$LinesOnAnnotation = scalar@Annotation;

for ($i=1;$i<$LinesOnAnnotation;$i++){
   @AnnotationData = split(",",$Annotation[$i]);
   chomp@AnnotationData;
   $Feature            = $AnnotationData[0];
   chomp$Feature;
   $Feature =~ s/\r//g;
   $Feature =~ s/\s//g;
   
   chomp$AnnotationData[1];
   $AnnotationData[1] =~ s/\r//g;
   $AnnotationData[1] =~ s/\s//g;
   $Gene{$Feature}     = $AnnotationData[1];
   
   chomp$AnnotationData[2];
   $AnnotationData[2] =~ s/\r//g;
   $AnnotationData[2] =~ s/\s//g;
   $EC{$Feature}       = $AnnotationData[2];

   chomp$AnnotationData[3];
   $AnnotationData[3] =~ s/\r//g;
   $AnnotationData[3] =~ s/\s//g;
   $Function{$Feature} = $AnnotationData[3];
}

open (TEST, ">$TestReport");
for ($i=0;$i<$LinesOnTrainingFile;$i++){
      if($i==0){
         print TEST "ORF,Gene,EC_Number,Initial_Annotation,All(%),";
         foreach $Class(@Classes){
            print TEST "$Class(%),";
         }
         foreach $Class(@Classes){
            print TEST "$Class,";
         }
         print TEST "SignificantClass,StatFreq\n";
         
      }else{
         $Feature = $TrainingMatrix[$i][0];
         print TEST "$Feature, $Gene{$Feature},$EC{$Feature},$Function{$Feature},$PercentageOfHit{$Feature},";
         
         @Percentages = ();
         @ChiValues = ();
         for ($j=0;$j<$nClasses*2+1;$j++){ 
            if ($j<$nClasses+1){
               if ($j>0){
                  $Percentage = $Percentages -> [$i][$j];
                  push @Percentages, "$Percentage";
                  print TEST "$Percentage,";
               }
            }else{
               $ChiValue = $Report -> [$i][$j-$nClasses];
               push @ChiValues, "$ChiValue";
               print TEST "$ChiValue,";
            }
         }
         if ( any { $_ > $Threshold} @ChiValues){
            push @SignificantGenes, "$Feature";
            foreach $ChiValue(@ChiValues){
               push (@{$ChiValues{$Feature}}, $ChiValue);
            }
            foreach $Percentage(@Percentages){
               push (@{$Percentages{$Feature}}, $Percentage);
            }
         }
         #my @TempChiValues = @ChiValues;
         #my @TempPercentages = @Percentages;
         #for (1..3){
         #   $InformativeIndex = first_index{$_ eq max@TempChiValues} @TempChiValues;
         #   $InformativeClass{$Feature} = $Classes[$InformativeIndex];
         #   splice@TempChiValues,$InformativeIndex,1;
         #   if ($Percentages[$InformativeIndex] == max@TempPercentages){
         #      $StatCondition = "Mostly Present";
         #      push @Present, $Feature;
         #      $PresentClassCount{$InformativeClass{$Feature}} += 1;
         #   }elsif($Percentages[$InformativeIndex] == min@TempPercentages){
         #      $StatCondition = "Mostly Absent";
         #      push @Absent, $Feature;
         #      $AbsentClassCount{$InformativeClass{$Feature}} += 1;
         #   }else{
         #      $StatCondition = "Uncertain";
         #   }
         #   splice@TempPercentages,$InformativeIndex,1;
         #   print TEST "$InformativeClass{$Feature},$StatCondition";
         #   if ($_ < 3){
         #      print TEST ",";
         #   }
         #}
         #print TEST "\n";
         $InformativeIndex = first_index{$_ eq max@ChiValues} @ChiValues;
         $InformativeClass{$Feature} = $Classes[$InformativeIndex];
         if ($Percentages[$InformativeIndex] == max@Percentages){
            $StatCondition = "Mostly Present";
            push @Present, $Feature;
            $PresentClassCount{$InformativeClass{$Feature}} += 1;
         }else{
            $StatCondition = "Mostly Absent";
            push @Absent, $Feature;
            $AbsentClassCount{$InformativeClass{$Feature}} += 1;
         }
         print TEST "$InformativeClass{$Feature},$StatCondition\n";
      }
}
close TEST;
$nPresent = scalar@Present;
$nAbsent = scalar@Absent;
$nSignificant = scalar@SignificantGenes;

#Loading TestReport file
($LinesOnTestReport, $ColumnsOnTestReport, @TestReportMatrix) = Matrix($TestReport);
$M = $ColumnsOnTestReport-1;      
$TestReportHeader = join ',', @{$TestReportMatrix[0]} [0..$M];
for ($i=1; $i<$LinesOnTestReport; $i++){
   $Feature = $TestReportMatrix[$i][0];
   $TestReportData{$Feature} = join ',', @{$TestReportMatrix[$i]} [1..$M];
}

open ($BOOLEAN, ">$BooleanAccessory");
open ($SIGNIFICANTPRESENT, ">$BooleanMostlyPresent");
open ($SIGNIFICANTABSENT, ">$BooleanMostlyAbsent");

for $FILE ($BOOLEAN, $SIGNIFICANTPRESENT,$SIGNIFICANTABSENT){
   print $FILE $TrainingHeader, "\n";
   close $FILE;
}
    
open (SIGNIFICANT, ">$SignificantReport");
   print SIGNIFICANT $TestReportHeader, "\n";
close SIGNIFICANT;
   
foreach $SignificantGene(@SignificantGenes){ 
   open (BOOLEAN, ">>$BooleanAccessory");
      print BOOLEAN "$SignificantGene,$BooleanData{$SignificantGene}\n";
   close BOOLEAN;
   
   if (any{ $_ eq $SignificantGene} @Present){
      push @SignificantPresent, $SignificantGene;
      $SignificantPresentClassCount{$InformativeClass{$SignificantGene}} += 1;
      open (SIGNIFICANTPRESENT, ">>$BooleanMostlyPresent");
         print SIGNIFICANTPRESENT "$SignificantGene,$BooleanData{$SignificantGene}\n";
      close SIGNIFICANTPRESENT;
   }elsif(any{ $_ eq $SignificantGene} @Absent){
      push @SignificantAbsent, $SignificantGene;
      $SignificantAbsentClassCount{$InformativeClass{$SignificantGene}} += 1;
      open (SIGNIFICANTABSENT, ">>$BooleanMostlyAbsent");
         print SIGNIFICANTABSENT "$SignificantGene,$BooleanData{$SignificantGene}\n";
      close SIGNIFICANTABSENT;
   }
   
   open (SIGNIFICANT, ">>$SignificantReport");
      print SIGNIFICANT "$SignificantGene,$TestReportData{$SignificantGene}\n";
   close SIGNIFICANT;
}
        
# Building dot plot
print "\nPlotting Chi values...";
chdir($OutPath);

if ($MappingFile){
   %Color = Mapping($MappingFile);
}else{
   foreach $Class(@Classes){
      $Color{$Class} = RGB;
   }
   if (scalar@Columns > 1){
      foreach $Class(@Classes2){
         $Color{$Class} = RGB;
      }
      foreach $Class(@Classes3){
         $Color{$Class} = RGB;
      }
   }
}
        
# Dot plot all clases
if ($DotPlot eq "AllClasses"){
   open(RSCRIPT, ">$PlotRScript");
      print RSCRIPT 'library(ggplot2)' . "\n";
      print RSCRIPT "df <- read.csv(\"$TestReport\")" . "\n";
      print RSCRIPT 'Colours <- c(';
      for ($i=0; $i<$nClasses-1;$i++){
         $Class = $Classes[$i];
         print RSCRIPT "$Class = \"$Color{$Class}\",";
      }
      print RSCRIPT "$Classes[$#Classes] = \"$Color{$Classes[$#Classes]}\"";
      print RSCRIPT ')' . "\n";
      print RSCRIPT 'ggplot(df, aes(ORF))';
      foreach $Class(@Classes){
         print RSCRIPT "+ geom_point(aes(y=$Class,color=\"$Class\"))";
      }
      print RSCRIPT "+ scale_colour_manual(values=Colours)";
      if($Method eq "X2"){
         @ChiConfidences = (0.5,0.6,0.7,0.8,0.9,0.95,0.995,0.999);
         #@ChiConfidences = (0.999);
         foreach $ChiConfidence(@ChiConfidences){
            print RSCRIPT "+ geom_hline(aes(yintercept = qchisq($ChiConfidence, df=$nClasses-1), linetype=\"$ChiConfidence\"))";
            print RSCRIPT "+ geom_text(aes(0, qchisq($ChiConfidence, df=$nClasses-1), label= round(qchisq($ChiConfidence, df=$nClasses-1),digits=2)), size=2, vjust=-0.25, hjust=0, nudge_x = 1)";
         }
      }
      print RSCRIPT "+ labs(x=\"Features\", y=\"Chi Values\", title= \"$Test Values of $nFeature Accessory Genes of
$ProjectName\'s Dataset\", color=\"$SelectedClass\", linetype=\"Confidence Interval\")";
      if($N > 100){
         print RSCRIPT '+ theme(axis.text.x = element_text(angle = 90, size=4, hjust = 1))' . "\n";
      }
      print RSCRIPT "\n";
      print RSCRIPT "ggsave(\"$Plot\")" . "\n";                   
   close RSCRIPT;
   
   system ("R CMD BATCH $PlotRScript");
   #system ("rm $PlotRScript $OutPath/*.Rout $OutPath/Rplots.pdf");
}
        
## Dot plot for class
#   if ($DotPlot eq "ForClass"){
#      foreach $Class(@Classes){
#         my $ClassPlotRScript = $OutPath ."/". $Class . "_DotPlotScript.R";
#         my $ClassPlot = $OutPath ."/". $Test ."_". $Class . "_DotPlot.pdf";
#         open(FILE, ">$ClassPlotRScript");
#            print FILE 'library(ggplot2)' . "\n";
#            print FILE "df <- read.csv(\"$TestReport\")" . "\n";
#            print FILE 'ggplot(df, aes(Feature))';
#            print FILE "+ geom_point(aes(y=$Class,color=\"$Class\"))";
#            if($Method eq "X2"){
#               @ChiConfidences = (0.9,0.95,0.975,0.99,0.999);
#               foreach $ChiConfidence(@ChiConfidences){
#                  print FILE "+ geom_hline(aes(yintercept = qchisq($ChiConfidence, df=$nClasses-1), linetype=\"$ChiConfidence\"))";
#               }
#            }
#            print FILE "+ labs(x=\"Genomic Feature\", y=\"Chi_Values\", title= \"$Test\", color=\"Class\", linetype=\"Confidence Interval\")";
#            if($N > 100){
#               print FILE '+ theme(axis.text.x = element_text(angle = 90, size=4, hjust = 1))' . "\n";
#            }
#            print FILE "\n";
#            print FILE "ggsave(\"$ClassPlot\")" . "\n";        
#         close FILE;
#         
#         system ("R CMD BATCH $ClassPlotRScript");
#         system ("rm -r $ClassPlotRScript $OutPath/*.Rout $OutPath/Rplots.pdf");
#                }
#        }
print "Done!\n";

#HeatMaps
my ($File, $SideColor, $Title, $FontSize, $Out);
my ($n);
my (@GenePresence);
my (%ClassCountOfSignificants);

open ($PRESENT, ">$PresentHeatMapRScript");
open ($ABSENT, ">$AbsentHeatMapRScript");
open($CORRELATED, ">$CorrelatedHeatMapRScript");

my ($Color, $Color1, $Color2, $Elements);
my $PlotData = 'freq';# <---Values to plot in the Heatmap. Could be stat or freq
    
for $FILE ($PRESENT, $ABSENT, $CORRELATED){
   print $FILE "\n";
   print $FILE "library(gplots)" . "\n";
   print $FILE "library(RColorBrewer)" . "\n";

   if($FILE eq $CORRELATED){
      $Color = "red";
      $Color1 = $Color;
      $Color2 = "oldlace";
      $n = $nStrains;
      $Elements = "Strains";
      $SideColor = 'ColSideColor=ClassesCol';
      $BooleanFile = $BooleanAccessory;
      $Title = "$ProjectName\'s Clustering
      Based on the Presence/Absence Profiles of $nSignificant
      Significant $SelectedClass Associated Genes ($Test,CI=$Threshold)";
      $Out = $HeatMapOut . "_Correlated.pdf";
      print $FILE "pdf(\"$Out\",height=11,width=11)" . "\n";
      
      print $FILE 'Strains <- c(';
      for ($i=0; $i<$nStrains-1;$i++){
         print $FILE "\"$Strains[$i]\",";
      }
      print $FILE "\"$Strains[$#Strains]\")" . "\n";

      foreach my $Class(@Classes){
         print $FILE "$Class <- c(";
         for ($i=0; $i<scalar@{$Classes{$Class}}-1; $i++){
            print $FILE "\"${$Classes{$Class}}[$i]\",";
         }
         print $FILE "\"${$Classes{$Class}}[$#{$Classes{$Class}}]\")" . "\n";
      } 
      if (scalar@Columns > 1){
         foreach my $NdClass (keys %Classes2){
            print $FILE "$NdClass <- c(";
            for ($i=0; $i<scalar@{$Classes2{$NdClass}}-1; $i++){
               print $FILE "\"${$Classes2{$NdClass}}[$i]\",";
            }
            print $FILE "\"${$Classes2{$NdClass}}[$#{$Classes2{$NdClass}}]\")" . "\n";
         }
         foreach my $RdClass(keys %Classes3){
            print $FILE "$RdClass <- c(";
            for ($i=0; $i<scalar@{$Classes3{$RdClass}}-1; $i++){
               print $FILE "\"${$Classes3{$RdClass}}[$i]\",";
            }
            print $FILE "\"${$Classes3{$RdClass}}[$#{$Classes3{$RdClass}}]\")" . "\n";
         }
      }
      
      print $FILE "df <- read.csv(\"$BooleanFile\")" . "\n";
      print $FILE 'rnames <- df[,1]' . "\n";
      print $FILE "Matrix <- data.matrix(df[,2:ncol(df)])" . "\n";
      print $FILE "rownames(Matrix) <- rnames" . "\n";
      #print $FILE "DisMatrix <- 1-cor(Matrix)" . "\n";
      print $FILE "DisMatrix <- data.matrix(dist(scale(t(Matrix), center=TRUE, scale=TRUE), method=\"euclidean\",upper=TRUE,diag=TRUE))" . "\n";
      $Matrix ="DisMatrix";
       
   }else{
      $Elements = "Features";
      $SideColor = 'RowSideColor=ClassesCol';
      
      if($FILE eq $ABSENT){
         $PermutationsFile = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "MostlyAbsent_Permutations.txt";
         if ($PlotData eq 'test'){
            $Color1 = "oldlace";
            $Color = "blue";
         }else{
            $Color1= "blue";
            $Color = "oldlace";
         }
         $n = scalar@SignificantAbsent;
         $Color2 = $Color;
         $BooleanFile = $BooleanAbsent;
         if ($PlotData eq "test"){
            $Title = "$Test(CI=$Threshold) values of $n Mostly Absent
            $SelectedClass Associated Genes on
            $ProjectName\'s Dataset";
         }else{
            $Title = "Frequency of $n Mostly Absent
$SelectedClass-Associated Genes
on $ProjectName\'s Dataset";
         }
         @GenePresence = @SignificantAbsent;
         %ClassCountOfSignificants = %SignificantAbsentClassCount;
         $Out = $HeatMapOut . "_Mostly_Absent.pdf";
         
      }elsif($FILE eq $PRESENT){
         $PermutationsFile = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "MostlyPresent_Permutations.txt";
         $Color1 = "oldlace";
         $Color = "red";
         $n = scalar@SignificantPresent;
         $Color2 = $Color;
         $BooleanFile = $BooleanPresent; 
         if ($PlotData eq "test"){
            $Title = "$Test(CI=$Threshold) values of $n Mostly Present
            $SelectedClass Associated Genes on
            $ProjectName\'s Dataset";
         }else{
            $Title = "Frequency of $n Mostly Present
$SelectedClass-Associated Genes
on $ProjectName\'s Dataset";
         }
         @GenePresence = @SignificantPresent;
         %ClassCountOfSignificants = %SignificantPresentClassCount;
         $Out = $HeatMapOut . "_Mostly_Present.pdf";
      }
      print $FILE "pdf(\"$Out\",height=9,width=9)" . "\n";
      $n = scalar@GenePresence;
      
      print $FILE 'Features <- c(';
      for ($i=0; $i<$n-1;$i++){
         print $FILE "\"$GenePresence[$i]\",";
      }
      print $FILE "\"$GenePresence[$#GenePresence]\")" . "\n";
      foreach $Class(@Classes){
         if($ClassCountOfSignificants{$Class}){
            print $FILE "$Class <- c(";
            $ClassCount{$Class} = ();
            for ($i=0; $i<$n; $i++){
               $Feature = $GenePresence[$i];   
               if($InformativeClass{$Feature} eq $Class){
                  $ClassCount{$Class} += 1;  
                  if($ClassCount{$Class} < $ClassCountOfSignificants{$Class}){
                     print $FILE "\"$Feature\",";
                  }else{
                     print $FILE "\"$Feature\")" . "\n";
                  }
               }
            }
         }
      }
      print $FILE "ORF <- c(";
      for ($i=0; $i<scalar@Classes-1; $i++){
         $Class = $Classes[$i];
         print $FILE "\"$Class\",";
      }
      print $FILE "\"$Classes[$#Classes]\")" . "\n";
      for ($i=0; $i<$n; $i++){
         $Feature = $GenePresence[$i];
         $Feature =~ s/\s//g;
         $Feature =~s/\b//g;
         print $FILE "$Feature <- c(";
         if ($PlotData eq 'stat'){
            for ($j=0; $j<scalar@{$ChiValues{$Feature}}-1; $j++){
               print $FILE "${$ChiValues{$Feature}}[$j],";
            }
            print $FILE "${$ChiValues{$Feature}}[$#{$ChiValues{$Feature}}])" . "\n";
         }elsif($PlotData eq 'freq'){
            for ($j=0; $j<scalar@{$Percentages{$Feature}}-1; $j++){
               print $FILE "${$Percentages{$Feature}}[$j],";
            }
            print $FILE "${$Percentages{$Feature}}[$#{$Percentages{$Feature}}])" . "\n";
         }  
      }
      
      print $FILE "Functions <- c(";
      for ($i=0; $i<$n-1; $i++){
         $Feature = $GenePresence[$i];
         print $FILE "\'$Gene{$Feature}\',";
      }
      print $FILE "\'$Gene{$GenePresence[$#GenePresence]}\')" . "\n";
      
      print $FILE "df <- data.frame(";
      for ($i=0; $i<$n-1; $i++){
         $Feature = $GenePresence[$i];
         print $FILE "$Feature,";
      }
      print $FILE "$GenePresence[$#GenePresence])" . "\n";
      
      print $FILE 'rnames <- ORF' . "\n";
      print $FILE "Matrix <- data.matrix(df)" . "\n";
      print $FILE "rownames(Matrix) <- rnames" . "\n";
      print $FILE "Matrix <- t(Matrix)" . "\n";
      print $FILE "hc <- hclust(as.dist(1-cor(t(Matrix))))" . "\n";
      $Matrix = "Matrix";
   }
   
   if($n < 50){
      $FontSize = 0.7;
   }elsif($n > 49 && $n < 80){
      $FontSize = 0.5;
   }elsif($n > 79 && $n < 150){
      $FontSize = 0.2;
   }elsif($n > 149){
      $FontSize = 0.1;
   }
   
   print $FILE "ClassesCol <- rep(\"black\",$n)" . "\n";
   
   foreach my $Class(keys %Classes){
      if ($Class){
         print $FILE "ClassesCol[$Elements %in% $Class] <- \"$Color{$Class}\"" . "\n";
      }
   }
   #print $FILE "Colors <- colorRampPalette(c(\"$Color1\", \"$Color2\"))($n)" . "\n";
   
   if ($FILE eq $CORRELATED){
      print $FILE "Colors <- colorRampPalette(c(\"$Color1\", \"$Color2\"))($n)" . "\n";
      if (scalar@Columns > 1){
         print $FILE "Classes2Col <- rep(\"black\",$n)" . "\n";
         foreach my $NdClass(keys %Classes2){
            if($NdClass){
               print $FILE "Classes2Col[$Elements %in% $NdClass] <- \"$Color{$NdClass}\"" . "\n";
            }
         }
         print $FILE "Classes3Col <- rep(\"black\",$n)" . "\n";
         foreach my $RdClass(keys %Classes3){
            if($RdClass){
               print $FILE "Classes3Col[$Elements %in% $RdClass] <- \"$Color{$RdClass}\"" . "\n";
            }
         }
      }
   }else{
      print $FILE "Colors <- colorRampPalette(c(\"$Color1\", \"$Color2\"))(1000)" . "\n";
   }

   print $FILE "HeatMap <- heatmap.2($Matrix," . "\n";
   print $FILE "$SideColor," ."\n";
   print $FILE "main = (\"$Title\")," . "\n";
   print $FILE "xlab = \"$SelectedClass\"," . "\n";
   print $FILE 'key.title = NA,' . "\n";
   if ($FILE eq $CORRELATED){
      print $FILE "key.xlab = \"Distance\"," . "\n";
   }else{
      if ($PlotData eq "stat"){
         print $FILE "key.xlab = \"$Test Values\"," . "\n";
      }else{
         print $FILE "key.xlab = \"Gene Frequency\"," . "\n";
      }
   }
   print $FILE 'trace = "none",' . "\n";                          
   print $FILE "tracecol = \"white\"," . "\n";
   print $FILE 'dendrogram = "both",' ."\n";
   
   if ($FILE eq $CORRELATED){
      if (scalar@Columns > 1){
         print $FILE "RowSideColor = Classes2Col," ."\n";
         print $FILE "colRow = Classes3Col," ."\n";
      }
      #print $FILE "distfun = function($Matrix) dist($Matrix, method = \"maximum\")," . "\n";
      print $FILE "ylab = \"Region Of Origin\"," . "\n";
      print $FILE 'srtCol= 90,' ."\n";
      print $FILE "cexCol=0.6," ."\n";
      print $FILE "cexRow=0.6," ."\n";
      print $FILE 'keysize=0.7,' ."\n";  
   }else{
      print $FILE "Rowv=as.dendrogram(hc)," . "\n";
      #print $FILE "ylab = \"Genomic Feature\"," . "\n";
      print $FILE 'srtCol= 20,' ."\n";
      #print $FILE 'adjCol= c(0.5,1),' ."\n";
      #print $FILE 'adjRow= c(0.25,0.5),' ."\n";
      print $FILE 'cexCol=0.75,' ."\n";
      print $FILE "cexRow=$FontSize," ."\n";
      print $FILE 'keysize=1,' ."\n";  
      print $FILE 'lhei=c(0.2,1),' ."\n";
      print $FILE 'labRow= Functions,' ."\n";
   }

   print $FILE 'col=Colors)' . "\n";
   
   if ($FILE ne $CORRELATED){
      print $FILE "Permutations <- rownames($Matrix)[HeatMap\$rowInd]" . "\n";
      print $FILE "write.csv(Permutations, file = \"$PermutationsFile\")" . "\n";
   }
   print $FILE 'HeatMap' . "\n";
                        
   print $FILE 'dev.off()';
   close $FILE;
}
system ("R CMD BATCH $CorrelatedHeatMapRScript");
system ("R CMD BATCH $PresentHeatMapRScript");
system ("R CMD BATCH $AbsentHeatMapRScript");
#system ("rm $OutPath/*.Rout $OutPath/*.R Rplots.pdf");
exit;