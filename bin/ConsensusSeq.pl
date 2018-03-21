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

my ($Usage, $ProjectName, $PresenceAbsence, $ConsensusSeq, $List, $CPUs,
    $MainPath, $ConsensusFile);

$Usage = "\tUsage: ConsensusSeq.pl <Main_Path> <Project Name> <ORFs Table> <Consensus Fasta File Name>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$MainPath      = $ARGV[0];
$ProjectName   = $ARGV[1];
$PresenceAbsence = $ARGV[2];
$ConsensusSeq = $ARGV[3];

my($Project, $ORFsPath, $TotalPresenceAbsence, $Row, $ORF, $ORFHmm, $LogFile,
   $LinesOnPresenceabsenceFile, $ColumnsOnPresenceabsenceFile);
my($i, $j);
my(@PresenceAbsence, @PresenceAbsenceFields, @PresenceAbsenceArray,
   @PresenceAbsenceMatrix);

$Project            = $MainPath ."/". $ProjectName;
$ORFsPath           = $Project ."/". "ORFs";

$LogFile            = $Project ."/". $ProjectName . ".log";

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

print "Loading the Presence/Absence file...";
($LinesOnPresenceabsenceFile, $ColumnsOnPresenceabsenceFile, @PresenceAbsenceMatrix) = Matrix($PresenceAbsence);
print "Done!\n";

print "Building a consensus multifasta file:\n";
for ($i=1; $i<$LinesOnPresenceabsenceFile; $i++){
        $ORF = $PresenceAbsenceMatrix[$i]->[0];
        $ORFHmm = $ORFsPath ."/". $ORF ."/". $ORF . ".hmm";
        
        system("hmmemit -c $ORFHmm >> $ConsensusSeq");
        
        Progress($LinesOnPresenceabsenceFile, $i);
}
exit;