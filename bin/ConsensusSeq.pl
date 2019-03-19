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

my ($Usage, $ProjectName, $Prefix, $PresenceAbsence, $ORFsPath, $OutPath);

$Usage = "\tUsage: ConsensusSeq.pl <Main_Path> <Project Name> <ORFs Table> <Consensus Fasta File Name>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName     = $ARGV[0];
$Prefix          = $ARGV[1];
$PresenceAbsence = $ARGV[2];
$ORFsPath        = $ARGV[3];
$OutPath         = $ARGV[4];

my($TotalPresenceAbsence, $Row, $ORF, $ORFHmm, $ConsensusSeq,
   $LinesOnPresenceabsenceFile, $ColumnsOnPresenceabsenceFile);
my($i, $j);
my(@PresenceAbsence, @PresenceAbsenceFields, @PresenceAbsenceArray,
   @PresenceAbsenceMatrix);

$ConsensusSeq = $OutPath ."/". $ProjectName . "_Consensus_" . $Prefix . ".fasta";

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