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

my ($Usage, $ProjectName, $AnnotationPath, $MainPath, $PresenceAbsence);

$Usage = "\tUsage: ORFsFunctions.pl <Project Name> <Annotation_Path> <Main_Path> <Presence-Absence>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName     = $ARGV[0];
$AnnotationPath  = $ARGV[1];
$MainPath        = $ARGV[2];
$PresenceAbsence = $ARGV[3];

my ($Project, $AnnotatedPresenceAnsence, $LinesOnPresenceAbsence,
    $ColumnsOnPresenceAbsence, $ORF, $cmd, $OTU, $AnnotationFile, $Function);
my ($i,$j);
my (@PresenceAbsenceMatrix, @Annotation);
my $Annotated = [ ];

$Project                  = $MainPath ."/". $ProjectName;
#$PresenceAbsence          = $Project ."/". $ProjectName . "_Presence_Absence.csv";
$AnnotatedPresenceAnsence = $Project ."/". $ProjectName . "_Annotated_Presence_Absence.csv";

($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, @PresenceAbsenceMatrix) = Matrix($PresenceAbsence);

$Annotated -> [0][0] = $PresenceAbsenceMatrix[0]->[0];

print "\nGetting annotation of each ORF:\n";
for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
   for ($j=1; $j<$ColumnsOnPresenceAbsence; $j++){
#for ($i=1; $i<10; $i++){
#   for ($j=1; $j<10; $j++){
      $OTU = $PresenceAbsenceMatrix[0][$j];
      $Annotated -> [$i][0] = $PresenceAbsenceMatrix[$i]->[0];
      $Annotated -> [0][$j] = $OTU;
      $ORF = $PresenceAbsenceMatrix[$i][$j];
      
      $AnnotationFile = $AnnotationPath ."/". $OTU ."/". $OTU . ".tsv";
     
      if ($ORF ne ""){
         #$cmd = `grep -r --include \"*.tsv\" \"$ORF\tCDS\" $AnnotationPath`;
         $cmd = `grep \"$ORF\tCDS\" $AnnotationFile`;
         
         @Annotation = split("\t",$cmd);
         $Function = $Annotation[$#Annotation];
         #$Function =~ s///g;
         chomp$Function;
         $Annotated -> [$i][$j] = $Function;
         
         print "\nThe function of the $PresenceAbsenceMatrix[$i]->[0] of $OTU is $Function";
      }else{
         $Annotated -> [$i][$j] = "";
      }
   }
   #Progress($LinesOnPresenceAbsence, $i);
}

print "\nBuilding Annotated file:\n";
open (FILE, ">$AnnotatedPresenceAnsence");
for ($i=0; $i<$LinesOnPresenceAbsence; $i++){
   for ($j=0; $j<$ColumnsOnPresenceAbsence; $j++){
#for ($i=0; $i<10; $i++){
   #for ($j=0; $j<10; $j++){
      print FILE $Annotated -> [$i][$j], ",";
   }
   print FILE "\n";
   #Progress($LinesOnPresenceAbsence, $i);
}
close FILE;
exit;


