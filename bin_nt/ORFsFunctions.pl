#!/usr/bin/perl -w
#################################################################################
#Scipt ContigsLength.pl                                                         #
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
    $ColumnsOnPresenceAbsence, $ORF, $cmd, $OTU, $AnnotationFile, $Function,
    $OTUORF, $Index, $ColumnsOnAnnotation, $Gene, $ECNumber, $ORFsFunctionsFile
    );
my ($i,$j);
my (@PresenceAbsenceMatrix, @Annotation, @PresenceAbsence, @Header, @ORFData);
my $Annotated = [ ];

$Project                  = $MainPath ."/". $ProjectName;
$AnnotatedPresenceAnsence = $Project ."/". $ProjectName . "_Annotated_Presence_Absence.csv";
$ORFsFunctionsFile        = $Project ."/". $ProjectName . "_ORFsAnnotation.csv";


print "\nGetting the function of ORFs:";

print "\nLoading the $PresenceAbsence file...";
@PresenceAbsence = ReadFile($PresenceAbsence);
$LinesOnPresenceAbsence = scalar@PresenceAbsence;

@Header = split(",",$PresenceAbsence[0]);
$ColumnsOnPresenceAbsence = scalar@Header;
print "Done!\n";

open (FILE, ">$ORFsFunctionsFile");
        print FILE "ORF,Gene,EC_Number,Product";
close FILE;

for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
        @ORFData = split (",",$PresenceAbsence[$i]);
        $ORF = $ORFData[0];
        chomp$ORF;
        shift@ORFData; 

        $Index = first_index { $_ ne "" } @ORFData;

        $OTUORF = $ORFData[$Index];
        $OTU = $Header[$Index+1];   
        
        $AnnotationFile = $AnnotationPath ."/". $OTU ."/". $OTU . ".tsv";

        $cmd = `grep \"$OTUORF\t\" $AnnotationFile`;
         
        @Annotation = split("\t",$cmd);
        $ColumnsOnAnnotation = $#Annotation;
        
        if(scalar@Annotation == 3){
            $Gene = "";
            $ECNumber = "";
        }elsif(scalar@Annotation == 4){
            if ($Annotation[2] =~ /^\d/){
                $ECNumber = $Annotation[2];
                chomp$ECNumber;
                $Gene = "";
            }else{
                $Gene     = $Annotation[2];
                chomp$Gene;
                $ECNumber = "";
            }
        }elsif(scalar@Annotation == 5){
            $Gene = $Annotation[2];
            $ECNumber = $Annotation[3];
        }
        
        $Function = $Annotation[$#Annotation];
        chomp$Function;
        $Function =~ s/,/-/g;
        
        open (FILE, ">>$ORFsFunctionsFile");
                print FILE "\n$ORF,$Gene,$ECNumber,$Function";
        close FILE;
        
        Progress($LinesOnPresenceAbsence, $i);
}

#($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, @PresenceAbsenceMatrix) = Matrix($PresenceAbsence);
#
#$Annotated -> [0][0] = $PresenceAbsenceMatrix[0]->[0];
#
#print "\nGetting annotation of each ORF:\n";
#for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
#   for ($j=1; $j<$ColumnsOnPresenceAbsence; $j++){
#      $OTU = $PresenceAbsenceMatrix[0][$j];
#      $Annotated -> [$i][0] = $PresenceAbsenceMatrix[$i]->[0];
#      $Annotated -> [0][$j] = $OTU;
#      $ORF = $PresenceAbsenceMatrix[$i][$j];
#      
#      $AnnotationFile = $AnnotationPath ."/". $OTU ."/". $OTU . ".tsv";
#     
#      if ($ORF ne ""){
#         #$cmd = `grep -r --include \"*.tsv\" \"$ORF\tCDS\" $AnnotationPath`;
#         $cmd = `grep \"$ORF\tCDS\" $AnnotationFile`;
#         
#         @Annotation = split("\t",$cmd);
#         $Function = $Annotation[$#Annotation];
#         $Function =~ s/,/-/g;
#         chomp$Function;
#         $Annotated -> [$i][$j] = $Function;
#         
#         #print "\nThe function of the $PresenceAbsenceMatrix[$i]->[0] of $OTU is $Function";
#      }else{
#         $Annotated -> [$i][$j] = "";
#      }
#   }
#   Progress($LinesOnPresenceAbsence, $i);
#}

#print "\nBuilding Annotated file:\n";
#open (FILE, ">$AnnotatedPresenceAnsence");
#for ($i=0; $i<$LinesOnPresenceAbsence; $i++){
#   for ($j=0; $j<$ColumnsOnPresenceAbsence; $j++){
#      print FILE $Annotated -> [$i][$j], ",";
#   }
#   print FILE "\n";
#   Progress($LinesOnPresenceAbsence, $i);
#}
#close FILE;

exit;


