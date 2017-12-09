#!/usr/bin/perl -w
#################################################################################
#Scipt ConsensusPanGenome.pl                                                    #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          19 de octubre de 2017                                           #
#################################################################################

use strict; 
use List::MoreUtils qw{any};
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $List, $CPUs, $MainPath);

$Usage = "\tUsage: ConsensusPanGenome.pl <Project Name> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$MainPath    = $ARGV[1];

my($Project, $ORFsPath, $PresenceAbsence, $TotalPresenceAbsence, $Row, $ORF,
   $ConsensusPanGenome, $ORFHmm, $LogFile);
my($i, $j);
my(@PresenceAbsence, @PresenceAbsenceFields, @PresenceAbsenceArray);

$Project = $MainPath ."/". $ProjectName;
$ORFsPath = $Project ."/". "ORFs";
$PresenceAbsence = $Project ."/". $ProjectName . "_Presence_Absence.csv";
$ConsensusPanGenome = $Project ."/". $ProjectName . "_Consensus_PanGenome.fasta";
$LogFile        = $Project ."/". $ProjectName . ".log";

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

@PresenceAbsence = ReadFile($PresenceAbsence);
$TotalPresenceAbsence = scalar@PresenceAbsence;

print "Loading the Presence/Absence file:\n";
for ($i=0; $i<$TotalPresenceAbsence; $i++){
     $Row = $PresenceAbsence[$i];
     @PresenceAbsenceFields = split(",",$Row);
     $j = scalar@PresenceAbsenceFields;
     push (@PresenceAbsenceArray, [@PresenceAbsenceFields]);
     Progress($TotalPresenceAbsence, $i);
}
print "Buiding a consensus Pan-Genome:\n";
for ($i=1; $i<$TotalPresenceAbsence; $i++){
        $ORF = $PresenceAbsenceArray[$i]->[0];
        $ORFHmm = $ORFsPath ."/". $ORF ."/". $ORF . ".hmm";
        
        system("hmmemit -c $ORFHmm >> $ConsensusPanGenome");
        
        Progress($TotalPresenceAbsence, $i);
}
exit;
