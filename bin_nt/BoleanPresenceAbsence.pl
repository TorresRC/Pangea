#!/usr/bin/perl -w
#################################################################################
#Scipt FormatPresenceAbsence.pl                                                 #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          11 de octubre de 2017                                           #
#################################################################################
use strict; 
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $List, $MainPath);

$Usage = "\tUsage: CoreGenome.pl <Project_Name> <List_File_Name> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List        = $ARGV[1];
$MainPath    = $ARGV[2];

my ($Project, $MainList, $PresenceAbsenceFile, $BoleanReport, $TotalQry,
    $LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, $Row, $Field,
    $BoleanInformativeReport, $Line, $nElements, $Count, $ORF);
my ($i, $j);
my (@List, @PresenceAbsence, @PresenceAbsenceFields, @PresenceAbsenceArray,
    @Elements, @PresenceAbsenceMatrix);
my $BoleanTable = [ ];
my $GenesAnnotationReport = [ ];

$Project                 = $MainPath ."/". $ProjectName;
$MainList                = $Project ."/". $List;
$PresenceAbsenceFile     = $Project ."/". $ProjectName . "_Presence_Absence.csv";
$BoleanReport            = $Project ."/". $ProjectName . "_Bolean_Presence_Absence.csv";
$BoleanInformativeReport = $Project ."/". $ProjectName . "_Bolean_Informative_Presence_Absence.csv";

print "\nLoading the Presence/Absence file:";
#@List = ReadFile($MainList);
#$TotalQry = scalar@List;
#@PresenceAbsence = ReadFile($PresenceAbsenceFile);
#$TotalPresenceAbsence = scalar@PresenceAbsence;
#
#for ($i=0; $i<$TotalPresenceAbsence; $i++){
#     $Row = $PresenceAbsence[$i];
#     @PresenceAbsenceFields = split(",",$Row);
#     chomp@PresenceAbsenceFields;
#     $o = scalar@PresenceAbsenceFields;
#     push (@PresenceAbsenceArray, [@PresenceAbsenceFields]);
#     Progress($TotalPresenceAbsence,$i);
#}

($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, @PresenceAbsenceMatrix) = Matrix($PresenceAbsenceFile);

print "Done!";

print "\nTraslating";

$BoleanTable -> [0][0] = $PresenceAbsenceMatrix[0][0];

for ($i=0; $i<$LinesOnPresenceAbsence; $i++){
     $BoleanTable -> [$i][0] = $PresenceAbsenceMatrix[$i][0];
     for ($j=0; $j<$ColumnsOnPresenceAbsence; $j++){
          $BoleanTable -> [0][$j] = $PresenceAbsenceMatrix[0][$j];  
     }
}

#Setting the bolean data ("1" for presence and "0" for absence)
for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
     for ($j=1; $j<$ColumnsOnPresenceAbsence; $j++){
          if (defined $PresenceAbsenceMatrix[$i]->[$j]){          
               $Field = $PresenceAbsenceMatrix[$i]->[$j];
               if ($Field ne ""){   
                    $BoleanTable -> [$i][$j] = "1";
               }else{
                    $BoleanTable -> [$i][$j] = "0";
               }
          }else{
               $BoleanTable -> [$i][$j] = "0";
          }
     }
     Progress($LinesOnPresenceAbsence,$i);
}

print "Writing file...";
open (FILE, ">$BoleanReport");
for ($i=0; $i<$LinesOnPresenceAbsence; $i++){
     for ($j=0; $j<$ColumnsOnPresenceAbsence; $j++){
          print FILE $BoleanTable -> [$i][$j], ",";
     }
     print FILE "\n";
     Progress($LinesOnPresenceAbsence,$i);
}
close FILE;

print "\nWriting a bolean file with only the informative ORFs...";
my @BoleanReport = ReadFile($BoleanReport);  
open (FILE,">$BoleanInformativeReport");
        print FILE "$BoleanReport[0]\n";
        for ($i=1; $i<scalar@BoleanReport;$i++){
                $Line = $BoleanReport[$i];
                @Elements = split(",",$Line);
                chomp@Elements;
                $ORF = $Elements[0];
                shift@Elements;
                $nElements = scalar@Elements;

                $Count = 0;
                foreach my $Element(@Elements){
                        if ($Element ne "0"){
                                $Count++;
                        }
                }
                if($Count < $nElements && $Count > 1){
                        print FILE "$Line\n";   
                }
        }
close FILE;

print "Done!\n";
exit;
