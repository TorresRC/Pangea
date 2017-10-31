#!/usr/bin/perl -w
#use strict;
use List::MoreUtils qw(uniq);
use lib "/home/rtorres/CoreGenome/src/lib";
use Routines;


my ($Usage, $ProjectName, $List, $CPUs, $MainPath);

$Usage = "\tUsage: ConsensusPanGenome.pl <Project Name> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
#$MainPath = $ARGV[1];

my($Project, $BoleanFileName, $MetaDataFileName, $ProbabilitiesFile, $nBoleanFile, $Line,
   $nBoleanFileFields, $N, $MetaData, $nMetaDataFile, $nMetaDataFileFields,
   $Region, $Strain, $Class, $nClasses, $Counter, $Hit, $Count, $Probe, $StrainHit,
   $StrainHits, $ProbeHit, $ProbeHits);
my($i, $j);
my(@BoleanFile, @BoleanFileFields, @BoleanTable, @MetaDataField, @MetaDataFile,
   @MetaDataFileFields, @MetaData, @Classes, @Strains);
my(%StrainClass, %Classes, %pClasses, %cpClasses, %ClassHits, %ProbeClass);
my $BoleanTable = [ ];
my $Probabilities = [ ];

$MainPath = "/home/rtorres/CoreGenome";
$Project = $MainPath ."/". $ProjectName;
$BoleanFileName = $Project ."/". $ProjectName . "_Bolean_Presence_Absence2.csv";
$MetaDataFileName = $Project ."/". 'Metadata.csv';
$ProbabilitiesFile = $Project ."/". $ProjectName . "_Probabilities.csv";
$LogFile        = $Project ."/". $ProjectName . ".log";

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

@BoleanFile = ReadFile($BoleanFileName);
$nBoleanFile = scalar@BoleanFile;
for ($i=0; $i<$nBoleanFile; $i++){
	$Line = $BoleanFile[$i];
	@BoleanFileFields = split(",",$Line);
	push (@BoleanTable, [@BoleanFileFields]);
}
$nBoleanFileFields = scalar@BoleanFileFields;
$N = $nBoleanFileFields-1;

@MetaDataFile = ReadFile($MetaDataFileName);
$nMetaDataFile = scalar@MetaDataFile;
for ($i=0; $i<$nMetaDataFile; $i++){
	$Line = $MetaDataFile[$i];
	@MetaDataFileFields = split(",",$Line);
	push (@MetaData, [@MetaDataFileFields]);
}
$nMetaDataFileFields = scalar@MetaDataFileFields;

for ($i=1;$i<$nMetaDataFile;$i++){
	$Class = $MetaData[$i]->[1];
	push @Classes, $Class;
}

@Classes = uniq(@Classes);
$nClasses = scalar@Classes;

for ($i=0;$i<$nClasses;$i++){
	$Counter = 0;
	for ($j=1;$j<$nMetaDataFile;$j++){
		$Strain = $MetaData[$j]->[0];
		$Class = $MetaData[$j]->[1];
		$StrainClass{$Strain} = $Class;
		if($Class eq $Classes[$i]){
			$Counter++;
		}
	}
	$Classes{$Classes[$i]} = $Counter;
	$pClasses{$Classes[$i]} = $Counter/$N;
	$cpClasses{$Classes[$i]} = 1-$Counter/$N;
}
foreach my $Class(@Classes){
print "\nThe Class $Class has $Classes{$Class} elements, and its probability is $pClasses{$Class} while the complemented probability is $cpClasses{$Class}";
}

$Count = 0;
for ($i=1; $i<$nBoleanFile; $i++){
	for ($j=1; $j<$nBoleanFileFields; $j++){
		$Hit = $BoleanTable[$i][$j];
		if ($Hit != 0){
			$Count++;
		}
	}
}
print "\nThe total of Hits into the bolean table are $Count\n";

foreach $Class(@Classes){
	$StrainHits = 0; 
	for ($i=1;$i<$nBoleanFileFields; $i++){
		$Strain = $BoleanTable[0][$i];
		if ($StrainClass{$Strain} eq $Class){
			for ($j=1;$j<$nBoleanFile;$j++){
				$StrainHit = $BoleanTable[$j][$i];
				$StrainHits += $StrainHit;
			}
		}
	}
	$ClassHits{$Class} = $StrainHits;		
}
foreach my $Class(@Classes){
	print "\nThe Class $Class has $ClassHits{$Class} hits";
}

foreach $Class(@Classes){
	for ($i=1;$i<$nBoleanFile;$i++){
		$Probe = $BoleanTable[$i][0];
		for ($j=1;$j<$nBoleanFileFields; $j++){
			$Strain = $BoleanTable[0][$j];
			if ($StrainClass{$Strain} eq $Class){
				$ProbeClassHit = $BoleanTable[$i][$j];
				$ProbeClass{$Probe}{$Class} += $ProbeClassHit;
			}
		}
	}
}

$Probabilities -> [0][0] = "Gene";
for ($i=1;$i<$nBoleanFile;$i++){
	$Probe = $BoleanTable[$i][0];
	$Probabilities -> [$i][0] = $Probe;
	for ($j=0; $j<scalar@Classes;$j++){
		$Class = $Classes[$j];
		$Probabilities -> [0][$j+1] = $Class;
		$Probabilities -> [$i][$j+1] = $ProbeClass{$Probe}{$Class}/$Classes{$Class}*100;
		print "\nThe probe $Probe has $ProbeClass{$Probe}{$Class} hits in class $Class";
	}
}

open (FILE, ">$ProbabilitiesFile");
for ($i=0; $i<$nBoleanFile; $i++){
     for ($j=0; $j<scalar@Classes+1; $j++){
          print FILE $Probabilities -> [$i][$j], ",";
     }
     print FILE "\n";
}
close FILE;

print "\n";
exit;

