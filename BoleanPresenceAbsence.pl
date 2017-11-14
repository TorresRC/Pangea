#!/usr/bin/perl -w
#################################################################################
#Scipt FormatPresenceAbsence.pl                                                 #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          11 de octubre de 2017                                           #
#################################################################################

use strict; 
use lib '/Users/rc/CoreGenome/src/lib';
use Routines;

my ($Usage, $ProjectName, $List, $MainPath);

$Usage = "\tUsage: CoreGenome.pl <Project_Name> <List_File_Name> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List = $ARGV[1];
$MainPath = $ARGV[2];

my ($Project, $MainList, $PresenceAbsenceFile, $BoleanReport, $TotalQry,
    $TotalPresenceAbsence, $Row, $Field, $BoleanInformativeReport, $Line,
    $nColumns, $Count);
my ($o, $i, $j);
my (@List, @PresenceAbsence, @PresenceAbsenceFields, @PresenceAbsenceArray, @Columns);
my $BoleanTable = [ ];
my $GenesAnnotationReport = [ ];

#$MainPath = "/Users/rc/CoreGenome/DiseaseMx";
$Project = $MainPath ."/". $ProjectName;
$MainList = $Project ."/". $List;
$PresenceAbsenceFile = $Project ."/". $ProjectName . "_Presence_Absence.csv";
$BoleanReport = $Project ."/". $ProjectName . "_Bolean_Presence_Absence.csv";
$BoleanInformativeReport = $Project ."/". $ProjectName . "_Bolean_Informative_Presence_Absence.csv";

#Loading the Presence/Absence genes file
@List = ReadFile($MainList);
$TotalQry = scalar@List;
@PresenceAbsence = ReadFile($PresenceAbsenceFile);
$TotalPresenceAbsence = scalar@PresenceAbsence;

for ($a=0; $a<$TotalPresenceAbsence; $a++){
     $Row = $PresenceAbsence[$a];
     @PresenceAbsenceFields = split(",",$Row);
     chomp@PresenceAbsenceFields;
     $o = scalar@PresenceAbsenceFields;
     push (@PresenceAbsenceArray, [@PresenceAbsenceFields]);
}

print "\nFormating the Absence\/Presence report to a bolean file...";

$BoleanTable -> [0][0] = $PresenceAbsenceArray[0][0];

for ($i=0; $i<$TotalPresenceAbsence; $i++){
     $BoleanTable -> [$i][0] = $PresenceAbsenceArray[$i][0];
     for ($j=0; $j<$o; $j++){
          $BoleanTable -> [0][$j] = $PresenceAbsenceArray[0][$j];  
     }
}

#Setting the bolean data ("1" for presence and "0" for absence)
for ($i=1; $i<$TotalPresenceAbsence; $i++){
     for ($j=1; $j<$o; $j++){
          if (defined $PresenceAbsenceArray[$i]->[$j]){          
               $Field = $PresenceAbsenceArray[$i]->[$j];
               if ($Field ne ""){   
                    $BoleanTable -> [$i][$j] = "1";
               }else{
                    $BoleanTable -> [$i][$j] = "0";
               }
          }else{
               $BoleanTable -> [$i][$j] = "0";
          }
     }
}

Writing the bolean file
open (FILE, ">$BoleanReport");
for ($i=0; $i<$TotalPresenceAbsence; $i++){
     for ($j=0; $j<$o; $j++){
          print FILE $BoleanTable -> [$i][$j], ",";
     }
     print FILE "\n";
}
close FILE;

my @BoleanReport = ReadFile($BoleanReport);  
open (FILE,">$BoleanInformativeReport");
        print FILE "$BoleanReport[0]\n";
        for ($i=0; $i<scalar@BoleanReport;$i++){
                $Line = $BoleanReport[$i];
                @Columns = split(",",$Line);
                chomp@Columns;
                $nColumns = scalar@Columns;

                $Count = 0;
                foreach my $Element(@Columns){
                        if ($Element ne "0"){
                                $Count++;
                        }
                }
                if($Count < $nColumns){
                        print FILE "$Line\n";   
                }
        }
close FILE;




print "Done!\n";
exit;
