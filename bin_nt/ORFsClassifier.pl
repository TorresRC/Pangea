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
$Correlation    = $ARGV[7];  # <- on/off
$Sort           = $ARGV[8];  # <- on/off
$Clusters       = $ARGV[9];  # <- on/off
$Dendrogram     = $ARGV[10]; # <- off/row/column/both
$Runs           = $ARGV[11];
$Threshold      = $ARGV[12];
$AnnotatedPresenceAbsence = $ARGV[13];

my($Test, $Run, $TestReport, $PercentagesReport, $Plot, $HeatMap, $PlotRScript,
   $LinesOnTrainingFile, $nFeature, $Line, $ColumnsOnTrainingFile, $N,
   $LinesOnMetaDataFile, $ColumnsOnMetaDataFile, $PossibleClass, $Column, $Class,
   $nClasses, $Element, $GlobalHits, $Hit, $Feature, $iClass, $a, $b, $c, $d,
   $nConfusion, $ChiConfidence, $Round, $HeatMapRScript, $Matrix, $CombinedInformative,
   $TestInformative, $PresenceInformative, $InformativeFeatures, $InformativeLines,
   $CombinedReport, $TestReportLine, $CombinedReportLine, $BoleanInformative,
   $TrainingHeader, $PresenceReportLine, $SummaryReport, $InformativeIndex,
   $InformativeClass, $SuperHeatMapRScript, $Percentage, $LinesOnAnnotation);
my($i, $j);
my(@TrainingFile, @TrainingFileFields, @TrainingMatrix, @MetaDataFile,
   @MetaDataFileFields, @MetaDataMatrix, @Classes, @Elements, @ChiConfidence,
   @ChiConfidences, @TestReport, @Combined, @TestReportData, @Presence,
   @PresenceData, @Annotation, @AnnotationData);
my(%ClassOfElement, %Elements, %pClass, %cpClass, %ClassHits,
   %HitsOfFeaturesInClass, %TotalFeatureHits, %Test,%PercentageOfFeatureInClass,
   %Gene, %EC, %Function);
my(%a, %b, %c, %d);
#my $Report = [ ];
#my $Percentages = [ ];
#my $Combined = [ ];

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

for ($i=1;$i<$LinesOnMetaDataFile;$i++){
	$Class = $MetaDataMatrix[$i]->[$Column];
        #$Class =~ s/\s//g;
	push @Classes, $Class;
}

@Classes = uniq(@Classes);
$nClasses = scalar@Classes;

foreach $Run (1 .. $Runs){
        
my $Report = [ ];
my $Percentages = [ ];
my $Combined = [ ];
        
        $TestReport          = $OutPath ."/". $Test . "_Values" . ".csv";
        $PercentagesReport   = $OutPath ."/". "PresencePercentages" . ".csv";
        $CombinedReport      = $OutPath ."/". $Test . "_ValuesAndPresencePercentages" . ".csv";
        
        $TestInformative     = $OutPath ."/". $Test . "_ValuesOfInformativeFeatures" . ".csv";
        $PresenceInformative = $OutPath ."/". $Test . "_PercentagesOfInformativeFeatures" . ".csv";
        $CombinedInformative = $OutPath ."/". $Test . "_ValuesAndPresencePercentagesOfInformativeFeatures" . ".csv";
        $BoleanInformative   = $OutPath ."/". $Test . "_BoleanInformativeFeatures" . ".csv";
        
        $SummaryReport       = $OutPath ."/". $Test . "_AssociatedGenes" . ".csv";
        
        $Plot                = $OutPath ."/". $Test . "_DotPlot" . ".pdf";
        $HeatMap             = $OutPath ."/". $Test . "_HeatMap" . ".png";
        $PlotRScript         = $OutPath ."/". "DotPlotScript" . ".R";
        $HeatMapRScript      = $OutPath ."/". "HeatMapScript" . ".R";
        $SuperHeatMapRScript = $OutPath ."/". "SuperHeatMapScript" . ".R";
        
        # Loading the bolean training file
        ($LinesOnTrainingFile, $ColumnsOnTrainingFile, @TrainingMatrix) = Matrix($TrainingFile);
        $nFeature = $LinesOnTrainingFile-1;
        $N = $ColumnsOnTrainingFile-1;
        
        $TrainingHeader = join ',', @{$TrainingMatrix[0]} [0..$N];
        #chop $TrainingHeader;
        #chop $TrainingHeader;
        
        open (FILE, ">$BoleanInformative");
                print FILE "$TrainingHeader\n";
        close FILE;
        
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
                        #print "\n$HitsOfFeaturesInClass{$Feature}{$Class})\t";
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
        
              $a= (($HitsOfFeaturesInClass{$Feature}{$Class}))+0.001; # hits de sonda a en clase a
              $b= (($TotalFeatureHits{$Feature}-$HitsOfFeaturesInClass{$Feature}{$Class}))+0.001; # Hits de sonda a que no están en clase A
              $c= (($Elements{$Class}-$HitsOfFeaturesInClass{$Feature}{$Class}))+0.001; # Numero de mismaches en clase A (numero de ceros en clase A)
              $d= ((($N-$Elements{$Class})-($TotalFeatureHits{$Feature}-$HitsOfFeaturesInClass{$Feature}{$Class})))+0.001; # Numero de ceros fuera de A
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
        
        ## Building output file
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
                $TestReportLine = $TestReport[$i];
                $PresenceReportLine = $Presence[$i];
                $CombinedReportLine = $Combined[$i];
                
                @PresenceData = split(",",$PresenceReportLine);
                shift@PresenceData;
                
                if ($i == 0){
                        open (CHI, ">$TestInformative");
                        open (PRESENCE, ">$PresenceInformative");
                        open (COMBINED, ">$CombinedInformative");
                        open (SUM, ">$SummaryReport");
                                print CHI $TestReportLine, "\n";
                                print PRESENCE $PresenceReportLine, "\n";
                                print COMBINED $CombinedReportLine, "\n";
                                print SUM "ORF,Gene,EC_Number,Product,";
                                foreach $Class(@Classes){
                           ##             #chop$Class;
                                        print SUM "$Class(%),";
                                }
                                print SUM "Represented_Class,Statistical_Prediction\n";
                        close CHI;
                        close PRESENCE;
                        close COMBINED;
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
                                        close CHI;
                                        close PRESENCE;
                                        close COMBINED;
                                        
                                        open (BOLEAN, ">>$BoleanInformative");
                                                my $InformativeFeature = `grep $Feature $TrainingFile`;
                                                chop $InformativeFeature;
                                                chop $InformativeFeature;
                                                print BOLEAN $InformativeFeature, "\n";
                                        close BOLEAN;
                                        
                                        $InformativeIndex = first_index{$_ eq max@TestReportData} @TestReportData;
                                        $InformativeClass = $Classes[$InformativeIndex];
                                        open (SUM, ">>$SummaryReport");
                                        $Function{$Feature} =~ s/\s//g;
                                        if ($PresenceData[$InformativeIndex] == max@PresenceData){
                                                print SUM "$Feature,$Gene{$Feature},$EC{$Feature},$Function{$Feature},";
                                                foreach $Percentage(@PresenceData){
                                                        print SUM "$Percentage,";
                                                }
                                                print SUM "$Classes[$InformativeIndex],Present\n";
                                        }else{
                                                #print SUM "$Feature,$Classes[$InformativeIndex],$PresenceData[$InformativeIndex],Absent\n";
                                                
                                                print SUM "$Feature,$Gene{$Feature},$EC{$Feature},$Function{$Feature},";
                                                foreach $Percentage(@PresenceData){
                                                        print SUM "$Percentage,";
                                                }
                                                print SUM "$Classes[$InformativeIndex],Absent\n";
                                        }
                                        close SUM;
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
                                                my $InformativeFeature = `grep $Feature $TrainingFile`;
                                                chop $InformativeFeature;
                                                chop $InformativeFeature;
                                                print BOLEAN $InformativeFeature, "\n";
                                        close BOLEAN;
                                }
                        }       
                }    
        }
        
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
                        
                        print RSCRIPT "png(\"$HeatMap\"," . "\n";
                        print RSCRIPT "width = 10*300," . "\n";
                        print RSCRIPT "height = 10*300," . "\n";
                        print RSCRIPT "res = 600," . "\n";
                        print RSCRIPT "pointsize = 5)" . "\n";
                        
                        print RSCRIPT 'Colors <- colorRampPalette(c("red", "yellow", "green"))(n=100)' . "\n";
                        
                        $Matrix = "Matrix";
                        print RSCRIPT "df <- read.csv(\"$TestInformative\")" . "\n";
                        print RSCRIPT 'rnames <- df[,1]' . "\n";
                        print RSCRIPT "$Matrix <- data.matrix(df[,2:ncol(df)])" . "\n";
                        print RSCRIPT "rownames($Matrix) <- rnames" . "\n";
                        
                        my $PresenceMatrix = "PresenceMatrix";
                        print RSCRIPT "TestV <- read.csv(\"$PresenceInformative\")" . "\n";
                        print RSCRIPT 'rnames <- TestV[,1]' . "\n";
                        print RSCRIPT "$PresenceMatrix <- data.matrix(TestV[,2:ncol(TestV)])" . "\n";
                        
                        if ($Correlation eq "on"){
                                print RSCRIPT "$Matrix <- cor($Matrix)" . "\n";
                        }
                        
                        print RSCRIPT "heatmap.2($Matrix," . "\n";
                        print RSCRIPT "main = \"$Test\"," . "\n";                        # Title
                        print RSCRIPT 'keysize = 0.5,' . "\n";
                        print RSCRIPT 'key.title = "Confidence",' . "\n";
                        print RSCRIPT 'key.xlab = "Significance",' . "\n";
                        #print RSCRIPT "cellnote = $PresenceMatrix," . "\n";
                        print RSCRIPT 'density.info="none",' . "\n";                     # Turns of density plot un legend
                        print RSCRIPT 'notecol = "black",' . "\n";                       # font of cell labels in black
                        print RSCRIPT 'trace = "none",' . "\n";                          # Turns of trace lines in heat map
                        
                        $Round = 5;
                        if ($Correlation eq "on"){
                                print RSCRIPT 'xlab = "Class",' . "\n";
                                print RSCRIPT 'ylab = "Class",' . "\n";
                                print RSCRIPT "cellnote = round($Matrix,$Round)," . "\n"; # Shows data in cell
                                print RSCRIPT 'srtRow= 90,' ."\n";
                                print RSCRIPT 'adjRow= c (0.5,1),' ."\n";
                                print RSCRIPT 'cexRow=0.8,' ."\n";
                        }else{
                                print RSCRIPT 'xlab = "Class",' . "\n";
                                print RSCRIPT 'ylab = "Feature",' . "\n";
                                if ($nClasses < 11 && $nFeature < 101){
                                        print RSCRIPT "cellnote = round($Matrix,$Round)," . "\n"; # Shows data in cell
                                }
                        }
                               
                        if ($Sort eq "off"){
                                print RSCRIPT 'Colv = "NA",' . "\n";                     # Turn off column sort
                                print RSCRIPT 'Rowv = "NA",' . "\n";                     # Turn off row sort
                        }elsif ($Sort eq "on"){
                                if ($Clusters eq "on"){
                                        #print RSCRIPT "distance = distfun($Matrix, method = \"manhattan\")," . "\n";
                                        #print RSCRIPT 'cluster = hclustfun(distance, method = "ward"),' . "\n";
                                        print RSCRIPT "Rowv = as.dendrogram(hclust(dist($Matrix, method = \"manhattan\"), method = \"ward\"))," . "\n";      # apply default clustering method just on rows
                                        print RSCRIPT "Colv = as.dendrogram(hclust(dist($Matrix, method = \"manhattan\"), method = \"ward\"))," . "\n";      # apply default clustering method just on columns
                                }
                        }
                        
                        if ($Dendrogram eq "off"){
                                print RSCRIPT 'dendrogram = "none",' ."\n";              # Hides dendrogram
                        }elsif($Dendrogram eq "row"){
                                print RSCRIPT 'dendrogram = "row",' ."\n";               # Shows dendrogram for rows
                        }elsif($Dendrogram eq "column"){
                                print RSCRIPT 'dendrogram = "column",' ."\n";            # Shows dendrogram for columns
                        }elsif($Dendrogram eq "both"){
                                print RSCRIPT 'dendrogram = "both",' ."\n";              # Shos dendrogram for rows and columns
                        }
                        
                        print RSCRIPT 'srtCol= 0,' ."\n";
                        print RSCRIPT 'adjCol= c (0.5,1),' ."\n";
                        print RSCRIPT 'cexCol=0.8,' ."\n";
                        print RSCRIPT 'col = Colors)' . "\n";                            # Use defined palette
                        
                        print RSCRIPT 'dev.off()';
                close RSCRIPT;
                system ("R CMD BATCH $HeatMapRScript");
                system ("rm $HeatMapRScript");
        }
           
        system ("rm $OutPath/*.Rout $OutPath/Rplots.pdf");
        
        print "Done!\n\n";
        
   
        # SuperHeat Map;
        #if ($HeatMapPlot eq "on"){
        #        open(RSCRIPT, ">$SuperHeatMapRScript");
        #                print RSCRIPT 'library(superheat)' . "\n";
        #                print RSCRIPT 'library(RColorBrewer)' . "\n";
        #                
        #                #print RSCRIPT "png(\"$HeatMap\"," . "\n";
        #                #print RSCRIPT "width = 10*300," . "\n";
        #                #print RSCRIPT "height = 10*300," . "\n";
        #                #print RSCRIPT "res = 300," . "\n";
        #                #print RSCRIPT "pointsize = 5)" . "\n";
        #                
        #                $Matrix = "Matrix";
        #                print RSCRIPT "df <- read.csv(\"$TestInformative\")" . "\n";
        #                print RSCRIPT 'rnames <- df[,1]' . "\n";
        #                print RSCRIPT "$Matrix <- data.matrix(df[,2:ncol(df)])" . "\n";
        #                print RSCRIPT "rownames($Matrix) <- rnames" . "\n";
        #                
        #                my $SummaryMatrix = "SummaryMatrix";
        #                print RSCRIPT "Summary <- read.csv(\"$SummaryReport\")" . "\n";
        #                print RSCRIPT 'rnames <- Summary[,1]' . "\n";
        #                print RSCRIPT "$SummaryMatrix <- data.matrix(Summary[,2:ncol(Summary)])" . "\n";
        #                print RSCRIPT "rownames($SummaryMatrix) <- rnames" . "\n";
        #                
        #                print RSCRIPT "superheat(as.matrix($Matrix)," . "\n";
        #                #print RSCRIPT "superheat($Matrix," . "\n";
        #                
        #                print RSCRIPT 'heat.pal = colorRampPalette(c("red", "yellow", "green"))(n=100),' . "\n";
        #                print RSCRIPT 'heat.na.col = "white",' . "\n";
        #                print RSCRIPT 'grid.vline.col = "white",' . "\n";
        #                
        #                #print RSCRIPT 'order.rows = "order(ORF)," . "\n";
        #        
        #                #print RSCRIPT 'yr = SummaryMatrix$Presence,' . "\n";
        #                #print RSCRIPT 'yr.plot.type = "bar",' . "\n";
        #                #print RSCRIPT 'yr.axis.name = "Percentage Of Presence on Associated Disease",' . "\n";
        #                #print RSCRIPT 'yr.plot.size = 0.2,' . "\n";
        #                #print RSCRIPT 'yr.point.size = 4,' . "\n";
        #                #print RSCRIPT 'yr.line.size = 2,' . "\n";
        #                
        #                print RSCRIPT 'left.label.size = 0.5,' . "\n";
        #                print RSCRIPT 'left.label.text.size = 3,' . "\n";
        #                print RSCRIPT 'bottom.label.size = 0.05,' . "\n";
        #                print RSCRIPT 'bottom.label.col = "white")' . "\n";
        #                #print RSCRIPT 'bottom.label.text.angle = 0,' . "\n";
        #                #print RSCRIPT 'bottom.label.text.alignment = "right")' . "\n";
        #                
        #                
        #                print RSCRIPT 'dev.off()';
        #        close RSCRIPT;
        #        system ("R CMD BATCH $SuperHeatMapRScript");
        #        system ("rm $SuperHeatMapRScript");
        #}
        
        $TrainingFile = $BoleanInformative;
        
}

exit;