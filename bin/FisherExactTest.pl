#!/usr/bin/perl -w

#################################################################################
#By:       Roberto C. Torres                                                    #
#e-mail:   torres.roberto.c@gmail.com                                           #
#################################################################################
use strict;
use FindBin;
use Statistics::R;
use lib "$FindBin::Bin/../lib";
use Routines;

my $MainPath = "$FindBin::Bin";
my $R = Statistics::R->new();

my ($PresenceAbsenceTable, $MetaDataFile, $QryGeneTable, $OutFile);

$QryGeneTable         = $ARGV[0];
$PresenceAbsenceTable = $ARGV[1];
$MetaDataFile         = $ARGV[2];
$OutFile              = $ARGV[3];

my ($LinesOnMetaDataFile, $ColumnsOnMetaDataFile, $Phenotype, $PhenotypeColumn,
    $PhenotypeName, $Strain, $Class, $ChosenClasses, $Control, $Case,
    $LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, $PresenceAbsenceHeader,
    $Gene, $LinesOnQry, $ColumnsOnQry, $pval, $ci, $or, $Pval, $CI, $OR, $CiInf,
    $CiSup, $CalcOR, $LogOdds, $Stdev, $WaldTest);
my ($i, $j, $nClasses, $nStrains, $nGenes, $N, $NumberOfStrainsWithGeneOnControl,
    $NumberOfStrainsWithGeneOnCase, $NumberOfStrainsWithoutGeneOnControl,
    $NumberOfStrainsWithoutGeneOnCase);
my (@MetaDataMatrix, @Strains, @Classes, @SelectedClasses,
    @PresenceAbsenceMatrix, @QryMatrix, @CI);
my (%ClassOfStrain, %StrainsInClass, %Classes, %BooleanData, %ClassHits,
    %HitsOfGeneInClass, %GenesInControl);

# Loading the metadata file
($LinesOnMetaDataFile, $ColumnsOnMetaDataFile, @MetaDataMatrix) = Matrix($MetaDataFile);


# Obtaining classes
print "\nThe following phenotypes were detected:";
for ($i=1;$i<$ColumnsOnMetaDataFile;$i++){
        $Phenotype = $MetaDataMatrix[0][$i];
        print "\n\t[$i] $Phenotype";
}
print "\n\nPlease type the index of the phenotype you want to analyze: ";
$PhenotypeColumn = <STDIN>;
chomp $PhenotypeColumn;

$PhenotypeName = $MetaDataMatrix[0]->[$PhenotypeColumn];

for ($i=1;$i<$LinesOnMetaDataFile;$i++){
   $Strain = $MetaDataMatrix[$i]->[0];
   push @Strains, $Strain;
   
   $Class = $MetaDataMatrix[$i]->[$PhenotypeColumn];
   $ClassOfStrain{$Strain} = $Class;
   $StrainsInClass{$Class}++; #   <--------------------------------------------- Number of elements in each class
   push @Classes, $Class;
   push (@{$Classes{$Class}}, $Strain);  
}

@Classes = uniq(@Classes);
$nClasses = scalar@Classes;
$nStrains = scalar@Strains;

print "\nThe following classes were detected in the phenotype $PhenotypeName:";
for ($i=0; $i<$nClasses; $i++){
    print "\n\t[$i] $Classes[$i]";
}
print "\n\nPlease select the two classes you want to compare:";
$ChosenClasses = <STDIN>;
chomp $ChosenClasses;

@SelectedClasses = split(",", $ChosenClasses);
$Control = $Classes[$SelectedClasses[0]];
$Case    = $Classes[$SelectedClasses[1]];

#Reading presence/absence table of all genes
($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, @PresenceAbsenceMatrix) = Matrix($PresenceAbsenceTable);
$nGenes = $LinesOnPresenceAbsence-1;
$N = $ColumnsOnPresenceAbsence-1;

$PresenceAbsenceHeader = join ',', @{$PresenceAbsenceMatrix[0]} [0..$N];

for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
   $Gene = $PresenceAbsenceMatrix[$i][0];
   $BooleanData{$Gene} = join ',', @{$PresenceAbsenceMatrix[$i]} [1..$N];
}


#Hits of each feature in each class
for ($i=1;$i<$LinesOnPresenceAbsence;$i++){
   $Gene = $PresenceAbsenceMatrix[$i][0];
   for ($j=1;$j<$ColumnsOnPresenceAbsence; $j++){
      $Strain = $PresenceAbsenceMatrix[0][$j];
      $HitsOfGeneInClass{$Gene}{$ClassOfStrain{$Strain}} += $PresenceAbsenceMatrix[$i][$j];# <- Total Feature Hits in Class
   } 
}

open (FILE, ">$OutFile");
print FILE "Gene,$Control(n=$StrainsInClass{$Control}),$Case(n=$StrainsInClass{$Case}),OR,LogOdds,Stdev,WaldTest,P-val,CI(95%)";
close FILE;

#Reading query table of all genes
($LinesOnQry, $ColumnsOnQry, @QryMatrix) = Matrix($QryGeneTable);

for ($i=1; $i<$LinesOnQry; $i++){
    $Gene = $QryMatrix[$i][0];
    
    $NumberOfStrainsWithGeneOnControl = $HitsOfGeneInClass{$Gene}{$Control};
    $NumberOfStrainsWithoutGeneOnControl = $StrainsInClass{$Control} - $HitsOfGeneInClass{$Gene}{$Control};
    $NumberOfStrainsWithGeneOnCase = $HitsOfGeneInClass{$Gene}{$Case};
    $NumberOfStrainsWithoutGeneOnCase = $StrainsInClass{$Case} - $HitsOfGeneInClass{$Gene}{$Case};
    
    #$CalcOR    = ($NumberOfStrainsWithGeneOnCase/$NumberOfStrainsWithGeneOnControl)/($NumberOfStrainsWithoutGeneOnCase/$NumberOfStrainsWithoutGeneOnControl);
    #$CalcLogOdds = log($CalcOR);
    
    $R->set('x', [$NumberOfStrainsWithGeneOnCase,$NumberOfStrainsWithoutGeneOnCase,$NumberOfStrainsWithGeneOnControl,$NumberOfStrainsWithoutGeneOnControl]);
    $R->run(
            q`fisher <- fisher.test(matrix(c(x), nrow=2))`,
            q`pval <- fisher["p.value"]`,
            q`ci <- fisher["conf.int"]`,
            q`or <- fisher["estimate"]`
            );
    $pval = $R -> get('pval');
    $ci   = $R -> get('ci');
    $or   = $R -> get('or');
    
    $R->stop();
    
    $Pval = ${$pval}[1];
    $CI = ${$ci}[0];
    $OR = ${$or}[3];
    
    @CI = split(" ",$CI);
    $CiInf = $CI[2];
    $CiSup = $CI[3];
    
    if ($NumberOfStrainsWithGeneOnCase > 0 && $NumberOfStrainsWithGeneOnControl > 0 && $NumberOfStrainsWithoutGeneOnCase > 0 && $NumberOfStrainsWithoutGeneOnControl){
      $Stdev = sqrt((1/$NumberOfStrainsWithGeneOnCase)+(1/$NumberOfStrainsWithGeneOnControl)+(1/$NumberOfStrainsWithoutGeneOnCase)+(1/$NumberOfStrainsWithoutGeneOnControl));
      $LogOdds = log($OR);
      $WaldTest  = $LogOdds/$Stdev;
    }else{
      $Stdev = "-";
      $LogOdds = "-";
      $WaldTest = "-";
    }
    
    open (FILE, ">>$OutFile");
    print FILE "\n$Gene,$NumberOfStrainsWithGeneOnControl,$NumberOfStrainsWithGeneOnCase,$OR,$LogOdds,$Stdev,$WaldTest,$Pval,$CiInf\-$CiSup";
    close FILE;
}

exit;










