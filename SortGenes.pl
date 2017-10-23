#!/usr/bin/perl -w
#################################################################################
#Scipt SortGenes.pl                                                             #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          14 de octubre de 2017                                           #
#################################################################################
use strict; 
use lib '/Users/rc/CoreGenome/src/lib';
use Routines;

my ($Usage, $ProjectName, $List, $TrustedORFeome, $MainPath);

$Usage = "\tUsage: SortGenes.pl <Project_Name> <List_File.ext> <Trusted_ORFeome.fasta> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List = $ARGV[1];
$TrustedORFeome = $ARGV[2];
$MainPath = $ARGV[3];

my ($MainList, $Project, $ORFeomesPath, $SeqExt, $LogFile, $ORFeome, $HeaderCharacter,
    $Gene, $GeneName, $GeneSeq, $GeneLength, $ContigsLengthReport, $Seq, $SortedORFeomePath,
    $SortedORFeome, $TrustedORFeomeFile, $SortedTrustedORFeome, $TrustedGene, $TrustedGeneName,
    $TrustedGeneSeq, $TrustedGeneLength);
my ($n, $i);
my (@List, @ORFeomes, @Genes, @TrustedGenes);
my (%GenesSeq, %GenesLength, %TrustedGenesLength, %TrustedGenesSeq);

#$MainPath = "/Users/rc/CoreGenome";
$Project = $MainPath ."/". $ProjectName;
$MainList = $Project ."/". $List;                                               
$ORFeomesPath = $MainPath ."/". "ORFeomes";
$TrustedORFeomeFile = $ORFeomesPath ."/". $TrustedORFeome;
$SortedORFeomePath = $ORFeomesPath ."/". "Sorted";
$SortedTrustedORFeome = $SortedORFeomePath ."/". $TrustedORFeome;
$SeqExt = ".ffn";
$LogFile = $Project ."/". $ProjectName . ".log";

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

MakeDir($SortedORFeomePath);
@ORFeomes = ReadFile($MainList);                                                     
$n = scalar@ORFeomes;


@TrustedGenes = ReadMultiFastaFile($TrustedORFeomeFile);
foreach $TrustedGene (@TrustedGenes){
        ($TrustedGeneName, $TrustedGeneSeq) = ReadSeq($TrustedGene);
        $TrustedGeneLength = length$TrustedGeneSeq;
        
        $TrustedGenesLength{"$TrustedGeneName"} = "$TrustedGeneLength";
        $TrustedGenesSeq{"$TrustedGeneName"} = "$TrustedGeneSeq";
}

open (FILE, ">$SortedTrustedORFeome");    
    foreach $TrustedGene (sort {$TrustedGenesLength{$b} <=> $TrustedGenesLength{$a} or $b cmp $a} keys %TrustedGenesLength){
            print FILE ">$TrustedGene\n$TrustedGenesSeq{$TrustedGene}\n";
    }
close FILE;


print "Sorting $n ORFeomes by Length: \n";

for ($i=0; $i<$n; $i++){
    $Seq = $ORFeomes[$i];
    %GenesSeq = ();
    %GenesLength = ();
    $ORFeome = $ORFeomesPath ."/". $Seq . $SeqExt;
    $SortedORFeome = $SortedORFeomePath ."/". $Seq . $SeqExt;
    
    @Genes = ReadMultiFastaFile($ORFeome);

    foreach $Gene (@Genes){
        ($GeneName, $GeneSeq) = ReadSeq($Gene);
        $GeneLength = length$GeneSeq;
        
        $GenesLength{"$GeneName"} = "$GeneLength";
        $GenesSeq{"$GeneName"} = "$GeneSeq";
    }

    open (FILE, ">$SortedORFeome");    
    foreach $Gene (sort {$GenesLength{$b} <=> $GenesLength{$a} or $b cmp $a} keys %GenesLength){
            print FILE ">$Gene\n$GenesSeq{$Gene}\n";
    }
    close FILE;
    Progress($n, $i);
}
exit;
    