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

#my $ConsensusScript = $MainPath ."/". "ConsensusSeq.pl";

my ($Usage, $TrainingFile, $MetadataFile, $OutPath, $Method, $X2, $IG, $OR,
    $PsCounts, $MI, $DotPlot, $HeatMapPlot, $Correlation, $Sort, $Clusters,
    $Dendrogram, $Runs, $Threshold, $AnnotatedPresenceAbsence, $ProjectName);

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
$ProjectName              = $ARGV[0];
$TrainingFile             = $ARGV[1];
$MetadataFile             = $ARGV[2];
$OutPath                  = $ARGV[3];
$AnnotatedPresenceAbsence = $ARGV[4];
$Method                   = $ARGV[5];  # <- X2, IG, OR, LO, MI
$DotPlot                  = $ARGV[6];  # <- AllClasses, ForClass
$Threshold                = $ARGV[7];

my($Test, $Run, $TestReport, $PercentagesReport, $Plot, $HeatMap, $PlotRScript,
   $LinesOnTrainingFile, $nFeature, $Line, $ColumnsOnTrainingFile, $N,
   $LinesOnMetaDataFile, $ColumnsOnMetaDataFile, $PossibleClass, $Column, $Class,
   $nClasses, $Element, $GlobalHits, $Hit, $Feature, $iClass, $a, $b, $c, $d,
   $nConfusion, $ChiConfidence, $Round, $HeatMapRScript, $Matrix, $CombinedInformative,
   $TestInformative, $PresenceInformative, $InformativeFeatures, $InformativeLines,
   $CombinedReport, $TestReportLine, $CombinedReportLine, $BoleanInformative,
   $TrainingHeader, $PresenceReportLine, $SignificantReport, $InformativeIndex,
   $InformativeClass, $SuperHeatMapRScript, $Percentage, $LinesOnAnnotation,
   $SelectedClass, $PresentSummaryReport, $AbsentSummaryReport, $BoleanInformativePresent,
   $BoleanInformativeAbsent, $TestInformativePresent, $TestInformativeAbsent,
   $StatCondition, $ChiValue, $SignificantGene, $BoleanSignificantGene, $SignificantGeneData,
   $Strain, $nStrains, $PresentHeatMap, $PresentHeatMapRScript, $AbsentHeatMap,
   $AbsentHeatMapRScript, $CorrelatedHeatMap, $CorrelatedHeatMapRScript, $BoleanFile,
   $HeatMapOut, $ConsensusSeq);
my($i, $j, $nPresent, $nAbsent, $BoleanPresent, $BoleanAbsent, $nSignificant);
my ($FILE, $PRESENT, $ABSENT, $CORRELATED, $BOLEAN);
my(@TrainingFile, @TrainingFileFields, @TrainingMatrix, @MetaDataFile,
   @MetaDataFileFields, @MetaDataMatrix, @Classes, @Elements, @ChiConfidence,
   @ChiConfidences, @TestReport, @Combined, @TestReportData, @Presence,
   @PresenceData, @Annotation, @AnnotationData, @Present, @Absent, @InformativeFiles,
   @Percentages, @ChiValues, @SignificantGenes, @Strains, @SignificantAbsent, @SignificantPresent);
my(%ClassOfElement, %Elements, %pClass, %cpClass, %ClassHits,
   %HitsOfFeaturesInClass, %TotalFeatureHits, %Test,%PercentageOfFeatureInClass,
   %Gene, %EC, %Function, %InformativeClass, %Color, %PresentClassCount,
   %AbsentClassCount, %ClassCount, %PercentageOfHit, %Phenotypes, %Classes, %SignificantPresentClassCount,
   %SignificantAbsentClassCount, %ChiValues);
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
print "\n\nPlease type the index of the desired class: ";
$Column = <STDIN>;
chomp $Column;

$SelectedClass = $MetaDataMatrix[0]->[$Column];

for ($i=1;$i<$LinesOnMetaDataFile;$i++){
   $Strain = $MetaDataMatrix[$i]->[0];
	$Class = $MetaDataMatrix[$i]->[$Column];
   push @Classes, $Class;
   push @Strains, $Strain;
   push (@{$Classes{$Class}}, $Strain);
}
@Classes = uniq(@Classes);
$nClasses = scalar@Classes;
$nStrains = scalar@Strains;
        
$TestReport               = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "ValuesOfAllGenes" . ".csv";
$SignificantReport        = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "ValuesOfSignificantGenes" . ".csv";
$BoleanInformative        = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "BoleanSignificantGenes" . ".csv";
$ConsensusSeq             = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "SignificantGenes" . ".fasta";
       
$Plot                     = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "ValuesOfAllGenes_DotPlot" . ".pdf";
$PlotRScript              = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "DotPlotScript" . ".R";
$PresentHeatMapRScript    = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "PresentHeatMapScript" . ".R";
$AbsentHeatMapRScript     = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "AbsentHeatMapScript" . ".R";
$CorrelatedHeatMapRScript = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI" ."_". "CorrelatedHeatMapScript" . ".R";
$HeatMapOut               = $OutPath ."/". $ProjectName ."_". $SelectedClass ."_". $Test . $Threshold . "CI";


@TrainingFile = ReadFile($TrainingFile);
        
# Loading the bolean training file
($LinesOnTrainingFile, $ColumnsOnTrainingFile, @TrainingMatrix) = Matrix($TrainingFile);
$nFeature = $LinesOnTrainingFile-1;
$N = $ColumnsOnTrainingFile-1;
        
$TrainingHeader = join ',', @{$TrainingMatrix[0]} [0..$N];
#chop $TrainingHeader;
#chop $TrainingHeader;
        
open (BOLEANINFORMATIVE, ">$BoleanInformative");
   print BOLEANINFORMATIVE "$TrainingHeader\n";
close BOLEANINFORMATIVE;
        
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

#Hits on each class
foreach $Class(@Classes){
   for ($i=1;$i<$ColumnsOnTrainingFile; $i++){
      $Element = "$TrainingMatrix[0][$i]";
      if ($ClassOfElement{$Element} eq $Class){
         for ($j=1;$j<$LinesOnTrainingFile;$j++){
            $ClassHits{$Class} += $TrainingMatrix[$j][$i];  # <- Total of hits in class
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
         $TotalFeatureHits{$Feature} += $TrainingMatrix[$i][$j]; # <- Total Feature Hits
         if ($ClassOfElement{$Element} eq $Class){
            $HitsOfFeaturesInClass{$Feature}{$Class} += $TrainingMatrix[$i][$j]; # <- Total Feature Hits in Class
         }
      }
      $PercentageOfHit{$Feature} = ($TotalFeatureHits{$Feature}/$N)*100;
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

print "\nBuilding OutPut Files...";
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
   $Function{$Feature} =~ s/\s//g;
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
         print TEST "Represented_Class,Statistical_Prediction\n";
         
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
         }        
         
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

@TestReport = ReadFile($TestReport);

open (BOLEAN, ">$BoleanInformative");
      print BOLEAN $TrainingFile[0], "\n";
close BOLEAN;

    
open (SIGNIFICANT, ">$SignificantReport");
print SIGNIFICANT $TestReport[0], "\n";
close SIGNIFICANT;
   
foreach $SignificantGene(@SignificantGenes){
   $BoleanSignificantGene = `grep -w $SignificantGene $TrainingFile`;
   $SignificantGeneData = `grep -w $SignificantGene $TestReport`;   
   open (BOLEAN, ">>$BoleanInformative");
      print BOLEAN $BoleanSignificantGene;
   close BOLEAN;
   
   if (any{ $_ eq $SignificantGene} @Present){
      push @SignificantPresent, $SignificantGene;
      $SignificantPresentClassCount{$InformativeClass{$SignificantGene}} += 1;
   }elsif(any{ $_ eq $SignificantGene} @Absent){
      push @SignificantAbsent, $SignificantGene;
      $SignificantAbsentClassCount{$InformativeClass{$SignificantGene}} += 1;
   }
   
   open (SIGNIFICANT, ">>$SignificantReport");
      print SIGNIFICANT $SignificantGeneData;
   close SIGNIFICANT;
}
print "Done!";

#system("perl $ConsensusScript $BoleanInformative $ConsensusSeq");
        
# Building dot plot
print "\nPlotting Chi values...";
chdir($OutPath);
        
# Dot plot all clases
if ($DotPlot eq "AllClasses"){
   open(RSCRIPT, ">$PlotRScript");
      print RSCRIPT 'library(ggplot2)' . "\n";
      print RSCRIPT "df <- read.csv(\"$TestReport\")" . "\n";
      print RSCRIPT 'ggplot(df, aes(ORF))';
      foreach $Class(@Classes){
         print RSCRIPT "+ geom_point(aes(y=$Class,color=\"$Class\"))";
      }
      if($Method eq "X2"){
         @ChiConfidences = (0.9,0.95,0.995,0.999);
         foreach $ChiConfidence(@ChiConfidences){
            print RSCRIPT "+ geom_hline(aes(yintercept = qchisq($ChiConfidence, df=$nClasses-1), linetype=\"$ChiConfidence\"))";
            print RSCRIPT "+ geom_text(aes(0, qchisq($ChiConfidence, df=$nClasses-1), label= round(qchisq($ChiConfidence, df=$nClasses-1),digits=2)), size=2, vjust=-0.25, hjust=0, nudge_x = 1)";
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
   #system ("rm $PlotRScript $OutPath/*.Rout $OutPath/Rplots.pdf");
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
         system ("rm -r $ClassPlotRScript $OutPath/*.Rout $OutPath/Rplots.pdf");
                }
        }
print "Done!\n";

#HeatMaps

my ($File, $SideColor, $Title, $FontSize, $Out);
my ($n);
my (@GenePresence);
my (%ClassCountOfSignificants);

open ($PRESENT, ">$PresentHeatMapRScript");
open ($ABSENT, ">$AbsentHeatMapRScript");
open($CORRELATED, ">$CorrelatedHeatMapRScript");

#my $t = scalar@{$Phenotypes{$Phenotype}};
#Random Color
#my @colors = map { 
#  join "", map { sprintf "%02x", rand(255) } (0..2) 
#} (0..63);
#
#$Color{$Phenotype} =

#################################################################################
$Color{None_Atrophic_Gastritis} = "darkgreen";
$Color{Atrophic_Gastritis} = "gold";
$Color{Intestinal_Metaplasia} = "orange";
$Color{Gastric_Cancer} = "red";
$Color{LatinAmerica} = "darkgreen";
$Color{Mexico} = "darkgreen";
#################################################################################

#################################################################################
my ($Country, $Population, $nCountries, $nPopulations, $Color1, $Color2, $Elements,
    $Color);
my (@Countries, @Populations);
my (%Countries, %Populations);
for ($i=1;$i<$LinesOnMetaDataFile;$i++){
   $Strain = $MetaDataMatrix[$i][0];
   $Country = $MetaDataMatrix[$i][3];
   $Population = $MetaDataMatrix[$i][5];
   push @Countries, $Country;
   push @Populations, $Population;
   push (@{$Countries{$Country}}, $Strain);
   push (@{$Populations{$Population}}, $Strain);
}
@Countries = uniq(@Countries);
@Populations = uniq(@Populations);
$nCountries = scalar@Countries;
$nPopulations = scalar@Populations;
#################################################################################

#################################################################################
#for ($i=0; $i<$n; $i++){
#   $Feature = $GenePresence[$i];
#   print $FILE "$Feature <- c(";
#   for ($j=0; $j<scalar@{$ChiValues{$Feature}}-1; $j++){
#      print $FILE "${$ChiValues{$Feature}}[$j],";
#   }
#   print $FILE "${$ChiValues{$Feature}}[$#{$ChiValues{$Feature}}])" . "\n";
#}
#
#foreach my $Country(keys %Countries){
#print $CORRELATED "$Country <- c(";
#   foreach (@{$Countries{$Country}}){
#      print $CORRELATED "\"$_\",";
#   }
#   print $CORRELATED "\"\")" . "\n";
#}
#foreach my $Population(keys %Populations){
#print $CORRELATED "$Population <- c(";
#   foreach (@{$Populations{$Population}}){
#      print $CORRELATED "\"$_\",";
#   }
#   print $CORRELATED "\"\")" . "\n";
#}
#################################################################################
    
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
      $BoleanFile = $BoleanInformative;
      $Title = "$ProjectName\'s Phylogeny
      Based on the Presence/Absence of $nSignificant
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
      
      print $FILE "df <- read.csv(\"$BoleanFile\")" . "\n";
      print $FILE 'rnames <- df[,1]' . "\n";
      print $FILE "Matrix <- data.matrix(df[,2:ncol(df)])" . "\n";
      print $FILE "rownames(Matrix) <- rnames" . "\n";
      print $FILE "DisMatrix <- 1-cor(Matrix)" . "\n";
      $Matrix ="DisMatrix";
       
   }else{
      $Color1 = "oldlace";
      $Elements = "Features";
      $SideColor = 'RowSideColor=ClassesCol';
      
      if($FILE eq $ABSENT){
         $Color = "blue";
         $n = scalar@SignificantAbsent;
         $Color2 = $Color;
         $BoleanFile = $BoleanAbsent;
         $Title = "$Test(CI=$Threshold) values of $n Mostly Absent
         $SelectedClass Associated Genes on
         $ProjectName Dataset";
         @GenePresence = @SignificantAbsent;
         %ClassCountOfSignificants = %SignificantAbsentClassCount;
         $Out = $HeatMapOut . "_Mostly_Absent.pdf";

         
      }elsif($FILE eq $PRESENT){
         $Color = "red";
         $n = scalar@SignificantPresent;
         $Color2 = $Color;
         $BoleanFile = $BoleanPresent;
         $Title = "$Test(CI=$Threshold) values of $n Mostly Present
         $SelectedClass Associated Genes on
         $ProjectName Dataset";
         @GenePresence = @SignificantPresent;
         %ClassCountOfSignificants = %SignificantPresentClassCount;
         $Out = $HeatMapOut . "_Mostly_Present.pdf";
      }
      print $FILE "pdf(\"$Out\",height=9,width=9)" . "\n";
      
      $n = scalar@GenePresence;
      
      if($n < 50){
         $FontSize = 1;
      }elsif($n > 49 && $n < 100){
         $FontSize = 0.8;
      }elsif($n > 99 && $n < 150){
         $FontSize = 0.5;
      }elsif($n > 149){
         $FontSize = 0.25;
      }
      
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
         print $FILE "$Feature <- c(";
         for ($j=0; $j<scalar@{$ChiValues{$Feature}}-1; $j++){
            print $FILE "${$ChiValues{$Feature}}[$j],";
         }
         print $FILE "${$ChiValues{$Feature}}[$#{$ChiValues{$Feature}}])" . "\n";
      }
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
   print $FILE "ClassesCol <- rep(\"black\",$n)" . "\n";
   
   #foreach my $Class(@Classes){
   #      print $FILE "$Class <- c(";
   #      for ($i=0; $i<scalar@{$Classes{$Class}}-1; $i++){
   #         print $FILE "\"${$Classes{$Class}}[$i]\",";
   #      }
   #      print $FILE "\"${$Classes{$Class}}[$#{$Classes{$Class}}]\")" . "\n";
   #   }
   
   foreach my $Class(keys %Classes){
      if ($Class){
         print $FILE "ClassesCol[$Elements %in% $Class] <- \"$Color{$Class}\"" . "\n";
      }
   }
   print $FILE "Colors <- colorRampPalette(c(\"$Color1\", \"$Color2\"))($n)" . "\n";
   
######################
      #print $FILE "CountryCol <- rep(\"black\",$n)" . "\n";
      #foreach my $Country(keys %Countries){
      #   print $FILE "CountryCol[$Elements %in% $Country] <- \"$Color{$Country}\"" . "\n";
      #}
      #print $FILE "PopulationCol <- rep(\"black\",$n)" . "\n";
      #foreach my $Population(keys %Populations){
      #   print $FILE "ClassesCol[$Elements %in% $Population] <- \"$Color{$Population}\"" . "\n";
      #}
######################

   print $FILE "heatmap.2($Matrix," . "\n";
   print $FILE "$SideColor," ."\n";
   print $FILE "main = (\"$Title\")," . "\n";
   print $FILE "xlab = \"$SelectedClass\"," . "\n";
   print $FILE 'key.title = NA,' . "\n";
   print $FILE "key.xlab = \"$Test Values\"," . "\n";
   print $FILE 'trace = "none",' . "\n";                          
   print $FILE "tracecol = \"white\"," . "\n";
   print $FILE 'dendrogram = "both",' ."\n";
   
   if ($FILE eq $CORRELATED){
      print $FILE 'srtCol= 90,' ."\n";
      print $FILE 'cexCol=0.6,' ."\n";
      print $FILE 'cexRow=0.6,' ."\n";
      print $FILE 'keysize=0.7,' ."\n";  
   }else{
      print $FILE "Rowv=as.dendrogram(hc)," . "\n";
      print $FILE "ylab = \"Genomic Feature\"," . "\n";
      print $FILE 'srtCol= 0,' ."\n";
      print $FILE 'adjCol= c(0.5,1),' ."\n";
      print $FILE 'adjRow= c(0.25,0.5),' ."\n";
      print $FILE 'cexCol=1,' ."\n";
      print $FILE "cexRow=$FontSize," ."\n";
      print $FILE 'keysize=1,' ."\n";
      print $FILE 'lhei=c(0.2,1),' ."\n";
   }

   print $FILE 'col=Colors)' . "\n";
                        
   print $FILE 'dev.off()';
   close $FILE;
}
system ("R CMD BATCH $CorrelatedHeatMapRScript");
system ("R CMD BATCH $PresentHeatMapRScript");
system ("R CMD BATCH $AbsentHeatMapRScript");
#system ("rm $OutPath/*.Rout $OutPath/*.R Rplots.pdf");
exit;