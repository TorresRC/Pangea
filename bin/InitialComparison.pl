#!/usr/bin/perl -w
#################################################################################
#Scipt InitialComparison.pl                                                     #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          April 11th 2017                                                 #
#################################################################################

use strict; 
use List::MoreUtils qw{any first_index};
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $Query, $Reference, $MolType, $eValue, $PIdent, $CPUs,
    $OutPath);

$Usage = "\tUsage: PreCoreGenome.pl <Project Name> <Genomes_List_File.ext> <TrustedORFeome.fasta> <Molecule type> <e-value> <Perc Ident> <CPUs> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$Query          = $ARGV[0];
$Reference      = $ARGV[1];
$MolType        = $ARGV[2];
$eValue         = $ARGV[3];
$PIdent         = $ARGV[4];
$CPUs           = $ARGV[5];
$OutPath        = $ARGV[6];

my ($QryPrefix, $RefPrefix, $SeqExt, $QryExt, $BlastPath, $BlastComparisons,
    $SamplesFeaturesPath, $FeaturesPath, $QryFile, $QryDb, $RefDb, $BlastReport,
    $LogFile, $RefFile, $Counter, $FTR, $FTRDir, $QryID, $RefID, $ObservedPIdent,
    $Index, $CatFTR, $QryOut, $RefOut, $Report, $RemainingQryFTRs);
my ($i);
my (@Lis, @QryIDs, @RefIDs, @BlastReport, @ReportFields, @AnalizedFTRs);
my (%FTR);

$QryPrefix = Prefix($Query);
$RefPrefix = Prefix($Reference);
$SeqExt = ".fasta";


if ($MolType eq "nucl"){
	$QryExt = ".ffn";
}elsif($MolType eq "prot"){
	$QryExt = ".faa";
}
     
$BlastPath           = $OutPath   ."/". "DataFiles" ."/". "Blast";
$BlastComparisons    = $BlastPath ."/". "Comparisons";
$SamplesFeaturesPath = $OutPath   ."/". "DataFiles" ."/". "SamplesFeatures";
$FeaturesPath        = $OutPath   ."/". "DataFiles" ."/". "Features";
$QryFile             = $SamplesFeaturesPath   ."/". $QryPrefix . $QryExt;
$RefFile             = $SamplesFeaturesPath   ."/". $RefPrefix . $QryExt;
$QryDb               = $BlastPath ."/". "DBs" ."/". $QryPrefix . "_BlastDb";
$RefDb               = $BlastPath ."/". "DBs" ."/". $RefPrefix . "_BlastDb";
$BlastReport         = $BlastComparisons ."/". $QryPrefix ."_Vs_". $RefPrefix ."_BlastReport.txt";
$LogFile             = $OutPath ."/". "Errors.log";
$Report              = $FeaturesPath ."/". "FeaturesMap.csv";

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

MakeDir($BlastComparisons);
MakeDir($FeaturesPath);

@QryIDs = AnnotatedGenes($QryFile);
@RefIDs = AnnotatedGenes($RefFile);

for($i=0;$i<@RefIDs;$i++){
    $Counter = sprintf "%.4d", $i+1;
    $FTR = "FTR_". $Counter;
    $FTRDir = $FeaturesPath ."/". $FTR;
    $FTR{$RefIDs[$i]} = $FTR;
    MakeDir($FTRDir);
}

Blast($MolType, $QryFile, $RefDb, $BlastReport, $eValue, $PIdent, $CPUs);
@BlastReport = ReadFile($BlastReport);

for ($i=0; $i<@BlastReport; $i++){
    @ReportFields   = SplitTab($BlastReport[$i]);
	$QryID          = $ReportFields[0];
	$RefID          = $ReportFields[1];
    $ObservedPIdent = $ReportFields[9];
    if ($PIdent < $ObservedPIdent){                      
        #Removing shared query features from the entire collection
        if ( any { $_ eq $QryID} @QryIDs){  
            $Index = first_index{$_ eq $QryID} @QryIDs;
            splice@QryIDs,$Index,1;
        }
        #Removing shared reference features from the entire collection
        if ( any { $_ eq $RefID} @RefIDs){  
            $Index = first_index{$_ eq $RefID} @RefIDs;
            splice@RefIDs,$Index,1;
        }
        $CatFTR = $FeaturesPath ."/". $FTR{$RefID} ."/". $FTR{$RefID} . $SeqExt;
        $QryOut = $FeaturesPath ."/". $FTR{$RefID} ."/". $QryPrefix ."_". $QryID . $SeqExt;        
        $RefOut = $FeaturesPath ."/". $FTR{$RefID} ."/". $RefPrefix ."_". $RefID . $SeqExt;                   
        if (any {$_ eq $RefID} @AnalizedFTRs){
            
        }else{               
            push @AnalizedFTRs, $RefID;
            Extract($QryPrefix,$QryDb,$MolType,$QryID,$QryOut);
            system("cat $QryOut >> $CatFTR");
            system("echo $FTR{$RefID},$QryPrefix,$QryID | cat >> $Report");
            if (!-e $RefOut){
                Extract($RefPrefix,$RefDb,$MolType,$RefID,$RefOut);
                system("cat $RefOut >> $CatFTR");
                system("echo $FTR{$RefID},$RefPrefix,$RefID | cat >> $Report");
            }
        }
    }
}

#Reference specific features
for($i=0; $i<@RefIDs; $i++){
    $CatFTR = $FeaturesPath ."/". $FTR{$RefIDs[$i]} ."/". $FTR{$RefIDs[$i]} . $SeqExt;
    $RefOut = $FeaturesPath ."/". $FTR{$RefIDs[$i]} ."/". $RefPrefix ."-". $RefIDs[$i] . $SeqExt;
    if (! -e $RefOut){
        Extract($RefPrefix,$RefDb,$MolType,$RefIDs[$i],$RefOut);
        system("cat $RefOut >> $CatFTR");
        system("echo $FTR{$RefIDs[$i]},$RefPrefix,$RefIDs[$i] | cat >> $Report");
    }
}

#Query specific features
#Solving query duplicates and fragmented features
####################################
$RemainingQryFTRs = $SamplesFeaturesPath ."/". "Remaining_" . $QryPrefix . "_Features" . $SeqExt;
for($i=0;$i<@QryIDs;$i++){
        $QryID = $QryIDs[$i];
        $QryOut = $FeaturesPath ."/". $QryPrefix ."_". $QryID . $SeqExt;
        $BlastReport = $BlastComparisons ."/". $QryPrefix ."_". $QryID ."_BlastReport.temp";
        Extract($QryPrefix,$QryDb,$MolType,$QryID,$QryOut);
        Blast($MolType, $QryOut, $RefDb, $BlastReport, $eValue, $PIdent, $CPUs);
        @BlastReport = ReadFile($BlastReport);
        if(@BlastReport){ 
            if ( any { $_ eq $QryID} @QryIDs){
                $Index = first_index{$_ eq $QryID} @QryIDs;
                splice@QryIDs,$Index,1;
            }
        }else{
            system("cat $QryOut >> $RemainingQryFTRs");
        }
        system ("rm $QryOut $BlastReport");
}

##Summary Report
#$PanGenomeSize=`grep ">" $PanGenomeSeq | wc -l`;
#chomp $PanGenomeSize;
#open (FILE, ">$Summary");
#        print FILE "Project Name: $ProjectName\n";
#        print FILE "Included Genomes List File: $List\n";
#        print FILE "Used CPU's: $CPUs\n";
#        print FILE "e-value Threshold: $eValue\n";
#        print FILE "Percentage of Identity: $PIdent\n";
#        print FILE "Number of taxa: $TotalQry\n";        
#        print FILE "Trusted ORFeome: $TotalTrustedORFs\n";
#close FILE;
#
##Progress File
#open (FILE, ">$Progress");
#	print FILE "NumberOfNewStrains,AnalyzedStrain,PanGenome,CoreGenome,NewGenes\n";
#	print FILE "0,$TrustedORFeomePrefix,$TotalTrustedORFs,$TotalTrustedORFs,0\n";
#	print FILE "1,$Qry,$PanGenomeSize,$CoreGenomeSize,$TotalNewORFs\n";
#close FILE;

exit;
