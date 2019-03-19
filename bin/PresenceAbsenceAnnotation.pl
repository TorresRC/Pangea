#!/usr/bin/perl -w
#################################################################################
#Scipt ORFsFunction.pl                                                          #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          10 de enero de 2018                                             #
#################################################################################
use strict; 
use List::MoreUtils qw{any first_index};
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $Prefix, $AnnotationPath, $OutPath, $PresenceAbsence);

$Usage = "\tUsage: ORFsFunctions.pl <Main_Path> <Project Name> <Annotation_Path> <Presence-Absence>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName       = $ARGV[0];
$Prefix            = $ARGV[1];
$AnnotationPath    = $ARGV[2];
$PresenceAbsence   = $ARGV[3];
$OutPath           = $ARGV[4];

my ($Project, $AnnotatedPresenceAnsence, $LinesOnPresenceAbsence,
    $ColumnsOnPresenceAbsence, $ORF, $cmd, $OTU, $AnnotationFile, $Function,
    $OTUORF, $Index, $ColumnsOnAnnotation, $Gene, $ECNumber, $ORFsFunctionsFile,
    $Line, $LocusTag
    );
my ($i,$j);
my (@PresenceAbsenceMatrix, @Annotation, @PresenceAbsence, @Header, @ORFData,
    @InFileName, @AnnotationFile);
my (%Annotation, %RepresentantLocus, %RepresentantOTU, %LocusTag, %Function, %RepresentantFunction, %Gene, %ECN,
    %RepresentantGene, %RepresentantECN);
my $Annotated = [ ];

$AnnotatedPresenceAnsence = $OutPath ."/". $ProjectName ."_". $Prefix . "_Annotated_Presence_Absence.csv";
        
print "\nLoading Presence Absence File...";
($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, @PresenceAbsenceMatrix) = Matrix($PresenceAbsence);
print "Done!";      

$Annotated -> [0][0] = $PresenceAbsenceMatrix[0]->[0];

print "\nGetting annotation of each ORF:\n";
for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
   for ($j=1; $j<$ColumnsOnPresenceAbsence; $j++){
      $OTU = $PresenceAbsenceMatrix[0][$j];
      $OTU =~s/\r//g;
      $Annotated -> [$i][0] = $PresenceAbsenceMatrix[$i]->[0];
      $Annotated -> [0][$j] = $OTU;
      $ORF = $PresenceAbsenceMatrix[$i][$j];
      
      $AnnotationFile = $AnnotationPath ."/". $OTU ."/". $OTU . ".tsv";
     
      if ($ORF){
         $cmd = `grep \"$ORF\t" $AnnotationFile`;
         
         @Annotation = split("\t",$cmd);
         $Function = $Annotation[$#Annotation];
         $Function =~ s/,/-/g;
         $Function =~ s/\b//g;
         $Function =~ s/ /_/g;
         $Function =~ s/\n//g;
         $Function =~ s/\s/_/g;
         $Function =~ s/\//-/g;
         chomp$Function;
         $Annotated -> [$i][$j] = $Function;

      }else{
         $Annotated -> [$i][$j] = "";
      }
   }
   Progress($LinesOnPresenceAbsence, $i);
}

print "\nBuilding Annotated file:\n";
open (FILE, ">$AnnotatedPresenceAnsence");
for ($i=0; $i<$LinesOnPresenceAbsence; $i++){
   for ($j=0; $j<$ColumnsOnPresenceAbsence; $j++){
      print FILE $Annotated -> [$i][$j], ",";
   }
   print FILE "\n";
   Progress($LinesOnPresenceAbsence, $i);
}
close FILE;
exit;