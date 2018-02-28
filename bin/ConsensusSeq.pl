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

my ($Usage, $ProjectName, $List, $CPUs, $MainPath, $ConsensusFile, $PanGenome,
    $StatORFs);

$Usage = "\tUsage: ConsensusSeq.pl <Main_Path> <Project Name> <ORFs Table> <Consensus Fasta File Name>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$MainPath      = $ARGV[0];
$ProjectName   = $ARGV[1];
$PanGenome     = $ARGV[2];
$StatORFs      = $ARGV[3];

my($Project, $ORFsPath, $PresenceAbsence, $TotalPresenceAbsence, $Row, $ORF,
   $ConsensusSeq, $ORFHmm, $LogFile, $LinesOnPresenceabsenceFile,
   $ColumnsOnPresenceabsenceFile);
my($i, $j);
my(@PresenceAbsence, @PresenceAbsenceFields, @PresenceAbsenceArray,
   @PresenceAbsenceMatrix);

$Project            = $MainPath ."/". $ProjectName;
$ORFsPath           = $Project ."/". "ORFs";

$LogFile            = $Project ."/". $ProjectName . ".log";

if ($PanGenome == 1){
        $PresenceAbsence    = $Project ."/". $ProjectName . "_Presence_Absence.csv";
        $ConsensusSeq       = $Project ."/". $ProjectName ."_". "Consensus_PanGenome" . ".fasta";
}elsif ($StatORFs == 1){
        $PresenceAbsence    = $Project ."/". "Statistically_AssociatedGenes.csv";
        $ConsensusSeq       = $Project ."/". $ProjectName ."_". "Statistically_AssociatedGenes" . ".fasta";
}

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

print "Loading the Presence/Absence file...";
($LinesOnPresenceabsenceFile, $ColumnsOnPresenceabsenceFile, @PresenceAbsenceMatrix) = Matrix($PresenceAbsence);
print "Done!\n";

print "Buiding a consensus multifasta file:\n";
for ($i=1; $i<$LinesOnPresenceabsenceFile; $i++){
        $ORF = $PresenceAbsenceMatrix[$i]->[0];
        $ORFHmm = $ORFsPath ."/". $ORF ."/". $ORF . ".hmm";
        
        system("hmmemit -c $ORFHmm >> $ConsensusSeq");
        
        Progress($LinesOnPresenceabsenceFile, $i);
}
exit;