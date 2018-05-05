#!/usr/bin/perl -w
#################################################################################
#Scipt MergeAnnotation.pl                                                       #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          09 de abril de 2018                                             #
#################################################################################
use strict;
use List::Util qw(sum);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;


my ($Usage, $ReportTable, $SourceName, $Source, $OutPut);

$Usage = "\tUsage: MergeAnnotation.pl <Report_File> <Source_File> <Output_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ReportTable = $ARGV[0];
$SourceName  = $ARGV[1];
$Source      = $ARGV[2];
$OutPut      = $ARGV[3];

my ($SeqName, $Annotation, $ReportData);
my ($LinesOnReport, $ColumnsOnReport, $LinesOnSource, $ColumnsOnSource,
    $ColumnsOnNewReport);
my ($i, $j);
my (@ReportMatrix, @SourceMatrix);
my (%Annotation);
my $Report = [ ];

print "\nLoading files...";
($LinesOnReport, $ColumnsOnReport, @ReportMatrix) = Matrix($ReportTable);
($LinesOnSource, $ColumnsOnSource, @SourceMatrix) = Matrix($Source);

$ColumnsOnNewReport = $ColumnsOnReport+1;
print "Done!";

print "\nRecovering the new annotation...";
for ($i=0; $i<$LinesOnSource; $i++){
    $SeqName = $SourceMatrix[$i][0];
    chomp $SeqName;
    $Annotation = $SourceMatrix[$i][1];
    chomp $Annotation;
    $Annotation =~ s/ $//g;
    $Annotation =~ s/\,/-/g;
    $Annotation =~s/\W/_/g;

    if ($SourceMatrix[$i][0] ne ""){
        $Annotation{$SeqName} = $Annotation;
    }else{
        $Annotation{$SeqName} = "";
    }
}
print "Done!";

print "\nMerging Annotation...";
for ($i=0; $i<$LinesOnReport; $i++){
    for ($j=0; $j<$ColumnsOnReport; $j++){
        $ReportData = $ReportMatrix[$i][$j];
        if ($ReportMatrix[$i][$j] ne ""){
            chomp $ReportData;
            $ReportData =~ s/\s//g;
            $Report -> [$i][$j] = $ReportData;
        }else{
            $Report -> [$i][$j] = "";
        }
    }
}

$Report -> [0][$ColumnsOnReport] = $SourceName;

for ($i=1; $i<$LinesOnReport; $i++){
    $SeqName = $ReportMatrix[$i][0];
    chomp $SeqName;
    if ($ReportMatrix[$i][0] ne ""){
        $Report -> [$i][$ColumnsOnReport] = $Annotation{$SeqName};
    }else{
        $Report -> [$i][$ColumnsOnReport] = "";
    }
}
print "Done!";

print "\nBuilding the final report...";
open (FILE, ">$OutPut");
for ($i=0; $i<$LinesOnReport; $i++){
     for ($j=0; $j<$ColumnsOnNewReport; $j++){
        my $Field = $Report -> [$i][$j];
        if (defined $Field){
        }else{
            $Field = "";
        }
            
        if ($j < $ColumnsOnNewReport){
            print FILE $Field, ",";
        }else{
            print FILE $Field;
        }
     }
     print FILE "\n";
}
close FILE;
print "Done!";
print "\n\tFinished!\n";
exit;
