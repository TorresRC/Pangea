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

my ($Usage, $ProjectName, $List, $CPUs, $MainPath, $ConsensusFile);

$Usage = "\tUsage: ConsensusSeq.pl <Main_Path> <Project Name> <Consensus File Name>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$MainPath      = $ARGV[0];
$ProjectName   = $ARGV[1];
$ConsensusFile = $ARGV[2];

my($Project, $ORFsPath, $PresenceAbsence, $TotalPresenceAbsence, $Row, $ORF,
   $ConsensusSeq, $ORFHmm, $LogFile, $LinesOnPresenceabsenceFile,
   $ColumnsOnPresenceabsenceFile);
my($i, $j);
my(@PresenceAbsence, @PresenceAbsenceFields, @PresenceAbsenceArray,
   @PresenceAbsenceMatrix);

$Project            = $MainPath ."/". $ProjectName;
$ORFsPath           = $Project ."/". "ORFs";
#$PresenceAbsence    = $Project ."/". $ProjectName . "_Presence_Absence.csv";
$PresenceAbsence    = $Project ."/". "Chi_Square_AssociatedGenes.csv";
$ConsensusSeq       = $Project ."/". $ProjectName ."_". $ConsensusFile . ".fasta";
$LogFile            = $Project ."/". $ProjectName . ".log";

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

#@PresenceAbsence = ReadFile($PresenceAbsence);
#$TotalPresenceAbsence = scalar@PresenceAbsence;

print "Loading the Presence/Absence file:\n";

($LinesOnPresenceabsenceFile, $ColumnsOnPresenceabsenceFile, @PresenceAbsenceMatrix) = Matrix($PresenceAbsence);

print "Buiding a consensus multifasta file:\n";
for ($i=1; $i<$LinesOnPresenceabsenceFile; $i++){
        $ORF = $PresenceAbsenceMatrix[$i]->[0];
        $ORFHmm = $ORFsPath ."/". $ORF ."/". $ORF . ".hmm";
        
        system("hmmemit -c $ORFHmm >> $ConsensusSeq");
        
        Progress($LinesOnPresenceabsenceFile, $i);
}


exit;
