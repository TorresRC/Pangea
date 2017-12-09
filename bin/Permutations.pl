#!/usr/bin/perl -w
#################################################################################
#Scipt ContigsLength.pl                                                         #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          11 de abril de 2017                                             #
#################################################################################

use strict; 
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $CPUs, $MainPath);

$Usage = "\tUsage: Permutations.pl <Project Name> <CPUs> <Main Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$CPUs        = $ARGV[1];
$MainPath    = $ARGV[2];

my($Project, $ORFeomesPath, $GenomesPath, $ORFsPath,
   $PreCoreGenomeFile, $PermutationsFile, $LogFile, $SeqExt, $HmmExt, $n, $Row,
   $o, $Count, $c, $ORF, $Qry, $QryOrfID, $ORFpath, $Hmm, $TempFile, $QryOrfeome,
   $Line, $BestHit, $Entry, $Strand, $Start, $End, $Issue, $d, $Coord1, $Coord2,
   $QryGenome, $e);
my(@PreCoreGenome, @Data, @FileIntoArray, @HmmReport, @BestHitData,);
my $Permutations = [ ];

$Project = $MainPath ."/". $ProjectName;

$ORFeomesPath      = $Project ."/". "ORFeomes";
$GenomesPath       = $Project ."/". "Genomes";
$ORFsPath          = $Project."/". "ORFs"; #Shared ORFs
$PreCoreGenomeFile = $Project ."/". "Pre-Core-Genome.csv";
$PermutationsFile  = $Project ."/". "Permutations.csv";

$LogFile           = $Project ."/". $ProjectName . ".log";

$SeqExt = ".fasta";
$HmmExt = ".hmm";

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

@PreCoreGenome = ReadFile($PreCoreGenomeFile);
$n = scalar@PreCoreGenome;

for ($a=0; $a<$n; $a++){
  $Row = $PreCoreGenome[$a];
  @Data = split(",",$Row);
  $o = scalar@Data;
  push (@FileIntoArray, [@Data]);
}

for ($b=0; $b<$o-1; $b++){
        $Qry = $FileIntoArray[0]->[$b+1];
        $Permutations -> [$b][0] = $Qry;
}

for ($b=1; $b<$o; $b++){
        $Qry = $FileIntoArray[0]->[$b];  
        $Count = 1;
        for ($c=1; $c<$n; $c++){
              $ORF = $FileIntoArray[$c]->[0];
              $QryOrfID = $FileIntoArray[$c]->[$b];
              
              $ORFpath = $ORFsPath ."/". $ORF;
              $Hmm = $ORFpath ."/". $ORF . $HmmExt;
              $TempFile = $ORFpath ."/". "temp";
              
              $QryGenome = $GenomesPath ."/". $Qry . $SeqExt;
              
              system("nhmmer -E 1e-30 --cpu $CPUs --noali --dfamtblout $TempFile $Hmm $QryGenome");
              
              @HmmReport = ReadFile($TempFile);
              open (FILE, ">>$TempFile");
              foreach $Line(@HmmReport){
                     print FILE "$Line\n";
                     }
              close FILE;
              @HmmReport = ReadFile($TempFile);
              
              if (-z $TempFile){
                     print "The ORF $ORF is not present in $Qry\n";
                     system("rm $TempFile");
                     next;
              }else{
                     $BestHit = $HmmReport[0];
                     $BestHit =~ s/\s^//g;
                     $BestHit =~ s/\s+/,/g;
                     @BestHitData = split(",",$BestHit);
                     
                     $Entry = $BestHitData[0];
                     $Strand = $BestHitData[8];
                     $Coord1 = $BestHitData[9];
                     $Coord2 = $BestHitData[10];
                     
                     if ($Coord1 < $Coord2){
                            $Start = $Coord1;
                            $End = $Coord2;
                     }else{
                            $Start = $Coord2;
                            $End = $Coord1;
                     }
                     
                     $Issue = $Start .",". $Strand .",". $Count;
                     $Count++;
              }
              $Permutations -> [$b][$c] = $Issue;  
              system("rm $TempFile");
       }
}
open (FILE, ">>$PermutationsFile");
for($d=0; $d<$o; $d++){
    for ($e=0; $e<$n; $e++){
        print FILE $Permutations -> [$d][$e], "|";
    }
    print FILE "\n";
}
close FILE;

exit;
