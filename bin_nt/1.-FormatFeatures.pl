#!/usr/bin/perl -w

################################################################################
#   Programa Format Features                                                    #
#   Nota: Se debe ajustar la ruta de lib y de la variable $PathSeq a las que    #
#   realmente se tengan donde se instalción el programa.                        #
#                                                                               #
# Programador:   Alfonso Méndez Tenorio, Roberto C. Torres                      #
# Fecha:         10 de enero de 2017                                            #
################################################################################
use strict; 
use lib '/home/bioinformatica/CoreGenome/src/lib';
#use lib '/Users/Roberto/Documents/lib';
use rutinas;



my ($Usage, $ProjectName, $List);

$Usage = "\tUsage: FormatFeatures.pl <Project Name> <List File Name>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List =$ARGV[1];

my ($MainPath, $Project, $FeaturesPath, $OriFeaturesPath,
    $FtdFeaturesPath, $OriExt, $FtdExt, $MainList, $Prefix, $QryFeature,
    $FtdFeature, $Header, $NoLines);

my ($FeatureCount, $FlagPrint, $Line, $SizeData, $FeatureStart, $FeatureEnd,
    $FeatureType, $FeatureStrand, $Temp, $FeatureGeneCount, $FeatureLen,
    $NewGeneFlag, $Prev, $PrevFlag, $GeneType, $FeatureGene, $FeatureGeneIndex,
    $LogFile);
my (@Data,);
my (@List, @Lines);

# if the program is installed in other location this variable should be changed to
# the new location
#$MainPath = "/Users/Roberto/CoreGenome";
$MainPath = "/home/bioinformatica/CoreGenome";
$Project = $MainPath ."/". $ProjectName;

#This directory must be present in the program path
#$FeaturesPath = $Project ."/". 'Features';
$FeaturesPath = $MainPath ."/". 'Features';

$MainList = $Project ."/". $List;

#This directory must be present in the program path
$OriFeaturesPath = $FeaturesPath ."/". "tbl";
$OriExt = ".tbl";
#$OriFeaturesPath = $FeaturesPath ."/". "txt";
#$OriExt = ".txt";
$FtdFeaturesPath = $FeaturesPath ."/". "ptt";
MakeDir($FtdFeaturesPath);
$FtdExt = ".ptt";
$LogFile = $Project ."/". $ProjectName . ".log";


open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

@List = ReadFile($MainList);

foreach my $File (@List){
    #$Prefix = Prefix($File);    
    $QryFeature = $OriFeaturesPath ."/". $File . $OriExt;
    $FtdFeature = $FtdFeaturesPath ."/". $File . $FtdExt;

    unless (open (ARCH, $QryFeature))
        {
            print "Can not open the $QryFeature file\n";
            print "Finished\n";
            exit;
        } 
        $Header = <ARCH>;   #Reads the title of the sequence
        @Lines = <ARCH>;  #Reads the content of the file
        #chomp (@lines);
        $NoLines = @Lines;
   
        print "$File:\n";
        print "Feature file has $NoLines lines\n";
    close ARCH;

    $FeatureCount = 0;
    $FlagPrint = 0;
    $FeatureGeneCount = 0;
        
    open (ARCH, ">$FtdFeature");
        print ARCH "$Header";
        print ARCH "$NoLines\n";
        print ARCH "Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product\n";

        $NewGeneFlag = 2;
        $Prev = -1;
        foreach $Line (@Lines) {
            if ($Line =~ /^[0-9]/) {
                $FeatureCount++;
                @Data = split(' ', $Line);
                $SizeData = @Data;
                $FeatureStart = $Data[0];
                $FeatureEnd = $Data[1];
            
                if ($SizeData > 2) {
                    $FeatureType = $Data[2]
                } else {
                    $FeatureType = "empty";  #pseudogen? 
                }
        
                if ($FeatureEnd > $FeatureStart) {
                    $FeatureStrand = '+';
                } else {
                    $Temp = $FeatureStart;
                    $FeatureStart = $FeatureEnd;
                    $FeatureEnd = $Temp;
                    $FeatureStrand = '-';
                }
        
                if ($FeatureType eq 'gene') {
                    $FeatureGeneCount++;
                    $NewGeneFlag = 1; 
                    $FeatureLen = $FeatureEnd-$FeatureStart+1;
                } else {
                    $NewGeneFlag = 2; #we are reading a new feature record but is not a gene    
                }
            } else {
                if ($NewGeneFlag ne 0) {
                    $PrevFlag = $NewGeneFlag;
                }
                @Data = split(' ', $Line);
                $GeneType = $Data[0];
                
                if (($FeatureType eq 'gene') && ($GeneType = 'gene') && ($NewGeneFlag == 1)) {
                    ($FeatureGene = $Data[1]);
                    $FeatureGene =~ s/\//-/g;
                    $NewGeneFlag = 0;   #we are not reading a new feature record
                    $FlagPrint = 1;
                }
            }

            if (($FlagPrint == 1) && ($FeatureType eq 'gene')){
        
                $FeatureGeneIndex = sprintf "%.4d", $FeatureGeneCount;
                $FeatureGene = $FeatureGene ."_". $FeatureGeneIndex;
        
                print ARCH "$FeatureStart..$FeatureEnd\t$FeatureStrand\t$FeatureLen\t$FeatureGeneCount\t$FeatureGene\t$FeatureGene\t$FeatureGene\t$FeatureGene\t$FeatureGene\n";
                $FlagPrint = 0;
                if (($FeatureGeneCount-$Prev) ne 1) {
                    print "Missing $FeatureGeneCount features\n";
                }
                $Prev = $FeatureGeneCount;
            } 
        }

    close ARCH;

    print "$FeatureCount features found\n";
    print "$FeatureGeneCount genes found\n\n";

}
print "\tFinished\n";
exit;