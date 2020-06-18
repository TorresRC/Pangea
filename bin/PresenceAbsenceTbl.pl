#!/usr/bin/perl -w
#################################################################################
#Scipt PresenceAbsenceTbl.pl                                                    #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          June 18th 2020                                                  #
#################################################################################

use strict; 
use List::MoreUtils qw{any true uniq first_index};
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $List, $PCore, $OutPath);

$Usage = "\tUsage: PresenceAbsenceTbl.pl <Mapping File> <Out Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$List        = $ARGV[0];
$PCore       = $ARGV[1];
$OutPath     = $ARGV[2];

my ($MappingFile, $Line,$FTR, $Sample, $PresenceAbsence, $CoreTbl,
    $CoreSize, $FTRPrevalence);
my ($i, $j);
my (@MapFile, @FTR, @FTRs, @List, @PresenceAbsence, @Line);
my (%FTR, %FTRPrevalence);
my $Report = [ ];

$MappingFile = $OutPath ."/". "DataFiles" ."/". "Features" ."/". "FeaturesMap.csv";
$PresenceAbsence = $OutPath ."/". "Initial_PresenceAbsence.csv";
$CoreTbl = $OutPath ."/". "Core_PresenceAbsence.csv";

@MapFile = ReadFile($MappingFile);
foreach $Line (@MapFile){
    @FTR = split(',',$Line);
    $FTR{$FTR[0]}{$FTR[1]} = $FTR[2];
    push@FTRs,$FTR[0];
}
@FTRs = sort { $a cmp $b } uniq(@FTRs);

@List = ReadFile($List);
$Report -> [0][0] = "Feature";
$Report -> [0][1] = "Prevalence";
for ($i=0;$i<@FTRs;$i++){
    for ($j=0;$j<@List;$j++){
        $FTR = $FTRs[$i];
        $Sample = $List[$j];    
        $Report -> [$i+1][0] = $FTR;
        $Report -> [0][$j+2] = $Sample;
        if ($FTR{$FTR}{$Sample}){
            $Report -> [$i+1][$j+2] = $FTR{$FTR}{$Sample};
            $FTRPrevalence{$FTR}++;
        }else{
            $Report -> [$i+1][$j+2] = "";
        } 
    }
    $Report -> [$i+1][1] = $FTRPrevalence{$FTR}/@List*100;
}

open (FILE, ">$PresenceAbsence");
for ($i=0;$i<@FTRs+1;$i++){
    for ($j=0;$j<@List+2;$j++){
        if ($j<@List+1){
            print FILE $Report -> [$i][$j], ",";
        }else{
            print FILE $Report -> [$i][$j];
        }
    }
    print FILE "\n";
}
close FILE;

exit;




