#!/usr/bin/perl -w
#################################################################################
#Scipt PresenceAbsenceTbl.pl                                                    #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          June 20th 2020                                                  #
#################################################################################

use strict; 
use List::MoreUtils qw{any true uniq first_index};
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $List, $PCore, $OutPath, $OutFile);

$Usage = "\tUsage: PresenceAbsenceTbl.pl <Mapping File> <Out Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$PCore       = $ARGV[0];
$OutPath     = $ARGV[1];
$OutFile     = $ARGV[2];

my ($PresenceAbsenceTbl,$LinesOnInputMatrix,$ColumnsOnInputMatrix,$UnalignedFTR,
    $AlignedFTR);
my ($i);
my (@InputMatrix);

$PresenceAbsenceTbl = $OutPath ."/". "Initial_PresenceAbsence.csv";
($LinesOnInputMatrix, $ColumnsOnInputMatrix, @InputMatrix) = Matrix($PresenceAbsenceTbl);

for ($i=1;$i<$LinesOnInputMatrix;$i++){
    if($InputMatrix[$i][1] >= $PCore){
        $UnalignedFTR = $OutPath ."/DataFiles/Features/". $InputMatrix[$i][0] ."/". $InputMatrix[$i][0] .".fasta"; 
        $AlignedFTR = $OutPath ."/DataFiles/Features/". $InputMatrix[$i][0] ."/". $InputMatrix[$i][0] .".aln.fasta";
        system("echo muscle -in $UnalignedFTR -out $AlignedFTR -quiet | cat >> $OutFile");
    }
}

exit;