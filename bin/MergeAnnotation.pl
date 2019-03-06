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


my ($Usage, $ReportTable, $Source, $Clean, $OutPut);

$Usage = "\tUsage: MergeAnnotation.pl <Report_File> <Source_File> <Output_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ReportTable = $ARGV[0];
$Source      = $ARGV[1];
$OutPut      = $ARGV[2];
$Clean       = $ARGV[3];

my ($SeqName, $Annotation, $ReportData, $Header, $Column, $Id, $Value);
my ($LinesOnReport, $ColumnsOnReport, $LinesOnSource, $ColumnsOnSource,
    $ColumnsOnNewReport, $NewColumns, $NewColumn);
my ($c, $i, $j);
my (@ReportMatrix, @SourceMatrix, @Columns);
my (%Metadata);
my $Report = [ ];

print "\nLoading files...";
($LinesOnReport, $ColumnsOnReport, @ReportMatrix) = Matrix($ReportTable);
($LinesOnSource, $ColumnsOnSource, @SourceMatrix) = Matrix($Source);

print "Done!";


# Obtaining headers
print "\nThe following columns were detected:";
for ($i=1;$i<$ColumnsOnSource;$i++){
        $Header = $SourceMatrix[0][$i];
        print "\n\t[$i] $Header";
}
print "\n\nPlease type the index of the desired column (more than one columns can be selected. e.g. 1,2,3): ";
$Column = <STDIN>;
chomp $Column;

@Columns = split(",",$Column);
chomp@Columns;
$NewColumns = scalar@Columns;
$ColumnsOnNewReport = $ColumnsOnReport+$NewColumns;

print "\nGetting metadata...";
foreach $c (@Columns){
    for ($i=0; $i<$LinesOnSource; $i++){
        $Id = $SourceMatrix[$i][0];
        chomp $Id;
        $Value = $SourceMatrix[$i][$c];
        chomp $Value;
        
        if ($Clean == 1){
            $Value =~ s/ $//g;
            $Value =~ s/\,/-/g;
            $Value =~s/\W/_/g;
        }
        
        if ($SourceMatrix[$i][$c] ne ""){
            $Metadata{$Id}{$c} = $Value;
        }else{
            $Metadata{$Id}{$c} = "";
        }
    }  
}
print "Done!";


print "\nMerging Annotation...";
for ($i=0; $i<$LinesOnReport; $i++){
    for ($j=0; $j<$ColumnsOnReport; $j++){
        $ReportData = $ReportMatrix[$i][$j];
        if (defined $ReportMatrix[$i][$j]){
            if ($ReportMatrix[$i][$j] ne ""){
                chomp $ReportData;
                $ReportData =~ s/\s//g;
                $Report -> [$i][$j] = $ReportData;
            }else{
                $Report -> [$i][$j] = "";
            }
        }
    }
}

for ($i=0;$i<$NewColumns;$i++){
    $Header = $SourceMatrix[0][$Columns[$i]];
    $NewColumn = $ColumnsOnReport+$i;
    $Report -> [0][$NewColumn] = $Header;
    for ($j=1; $j<$LinesOnReport; $j++){
        $Id = $ReportMatrix[$j][0];
        chomp $Id;
        $Report -> [$j][$NewColumn] = $Metadata{$Id}{$Columns[$i]};
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
