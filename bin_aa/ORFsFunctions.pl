#!/usr/bin/perl -w
#################################################################################
#Scipt ContigsLength.pl                                                         #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          10 de enero de 2018                                             #
#################################################################################

use strict; 
use List::MoreUtils qw{any};
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $AnnotationPath, $MainPath);

$Usage = "\tUsage: ORFsFunctions.pl <Project Name> <Annotation_Path> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName     = $ARGV[0];
$AnnotationPath  = $ARGV[1];
$MainPath        = $ARGV[2];

my ($Project, $PresenceAbsence, $AnnotatedPresenceAnsence, $LinesOnPresenceAbsence,
    $ColumnsOnPresenceAbsence, $ORF, $cmd);
my ($i,$j);
my (@PresenceAbsenceMatrix, @Annotation);
my $Annotated = [ ];

$Project                  = $MainPath ."/". $ProjectName;
$PresenceAbsence          = $Project ."/". $ProjectName . "_Presence_Absence_B.csv";
$AnnotatedPresenceAnsence = $Project ."/". $ProjectName . "_Annotated_Presence_Absence.csv";

($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, @PresenceAbsenceMatrix) = Matrix($PresenceAbsence);

$Annotated -> [0][0] = $PresenceAbsenceMatrix[0]->[0];

print "\nAnnotating:\n";
for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
   for ($j=1; $j<$ColumnsOnPresenceAbsence; $j++){
      $Annotated -> [$i][0] = $PresenceAbsenceMatrix[$i]->[0];
      $Annotated -> [$j][0] = $PresenceAbsenceMatrix[0][$j];
      $ORF = $PresenceAbsenceMatrix[$i][$j];
      
      if ($ORF ne ""){
         $cmd = `grep -r --include \"*.tsv\" \"$ORF\tCDS\" $AnnotationPath`;
         
         @Annotation = split("\t",$cmd);
         $Annotated -> [$i][$j] = $Annotation[3];
         
         print "\nThe function of the $PresenceAbsenceMatrix[$i]->[0] of $PresenceAbsenceMatrix[0][$j] is $cmd";
      }else{
         $Annotated -> [$i][$j] = "";
      }
   }
   #Progress($LinesOnPresenceAbsence, $i);
}

print "\nBuilding Annotated file:\n";
open (FILE, ">$AnnotatedPresenceAnsence");
for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
   for ($j=1; $j<$ColumnsOnPresenceAbsence; $j++){
      print FILE $Annotated -> [$i][$j], ",";
   }
   print FILE "\n";
   #Progress($LinesOnPresenceAbsence, $i);
}
close FILE;
exit;


