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


my ($PresenceAbsenceTable, $MetaDataFile, $QryGeneTable);


$QryGeneTable         = $ARGV[0];
$PresenceAbsenceTable = $ARGV[1];
$MetaDataFile         = $ARGV[2];


my ($LinesOnMetaDataFile, $ColumnsOnMetaDataFile, $Phenotype, $PhenotypeColumn,
    $PhenotypeName, $Strain, $Class, $ChosenClasses, $Control, $Case,
    $LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, $PresenceAbsenceHeader,
    $Gene, $LinesOnQry, $ColumnsOnQry);
my ($i, $j, $nClasses, $nStrains, $nGenes, $N, $NumberOfStrainsWithGeneOnControl,
    $NumberOfStrainsWithGeneOnCase, $NumberOfStrainsWithoutGeneOnControl,
    $NumberOfStrainsWithoutGeneOnCase);
my (@MetaDataMatrix, @Strains, @Classes, @SelectedClasses,
    @PresenceAbsenceMatrix, @QryMatrix);
my (%ClassOfStrain, %StrainsInClass, %Classes, %BooleanData, %ClassHits, %HitsOfGeneInClass,
    %GenesInControl);

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

#Reading query table of all genes
($LinesOnQry, $ColumnsOnQry, @QryMatrix) = Matrix($QryGeneTable);

for ($i=1; $i<$LinesOnQry; $i++){
    $Gene = $QryMatrix[$i][0];
    
    $NumberOfStrainsWithGeneOnControl = $HitsOfGeneInClass{$Gene}{$Control};
    $NumberOfStrainsWithoutGeneOnControl = $StrainsInClass{$Control} - $HitsOfGeneInClass{$Gene}{$Control};
    $NumberOfStrainsWithGeneOnCase = $HitsOfGeneInClass{$Gene}{$Case};
    $NumberOfStrainsWithoutGeneOnCase = $StrainsInClass{$Case} - $HitsOfGeneInClass{$Gene}{$Case};

    
print "\n$Gene". '<- matrix(c('."$NumberOfStrainsWithGeneOnCase,$NumberOfStrainsWithoutGeneOnCase,$NumberOfStrainsWithGeneOnControl,$NumberOfStrainsWithoutGeneOnControl)" . ',nrow=2)'."\n";

}








