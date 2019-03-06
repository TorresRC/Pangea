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

my ($Usage, $PresenceAbsence, $ORFsPath, $OutPath);

$Usage = "\tUsage: ConsensusSeq.pl <Main_Path> <Project Name> <ORFs Table> <Consensus Fasta File Name>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$PresenceAbsence = $ARGV[0];
$ORFsPath = $ARGV[1];
$OutPath = $ARGV[2];

my($TotalPresenceAbsence, $Row, $ORF, $ORFHmm, $LogFile, $ConsensusSeq,
   $LinesOnPresenceabsenceFile, $ColumnsOnPresenceabsenceFile);
my($i, $j);
my(@PresenceAbsence, @PresenceAbsenceFields, @PresenceAbsenceArray,
   @PresenceAbsenceMatrix);

#$LogFile            = $OutPath ."/". $ProjectName . ".log";

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

$ConsensusSeq = $OutPath ."/". "ConsensusSeq.fasta";

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