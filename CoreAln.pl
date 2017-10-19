#!/usr/bin/perl -w

#################################################################################
#   Programa CoreAln                                                       #
#   Nota: Se debe ajustar la ruta de lib y de la variable $PathSeq a las que    #
#   realmente se tengan donde se instalción el programa.                        #
#                                                                               #
# Programador:   Roberto C. Torres                                              #
# Fecha:         26 de junio de 2017                                            #
#################################################################################

use strict;
use lib '/home/bioinformatica/CoreGenome/src/lib';
use Routines;

my ($Usage, $ProjectName, $List, $CPUs);

$Usage = "\tUsage: CoreAln.pl <Project Name> <List> <CPUs>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List = $ARGV[1];
$CPUs = $ARGV[2];


my($MainPath, $Project, $MainList, $ORFsPath, $CoreGenomeFile, $CoreAlnPath,
   $CoreSeqPath, $LogFile, $SeqExt, $AlnExt, $HmmExt, $n, $m, $Row, $o, $ORF,
   $CoreORF, $AlnCoreORF, $TaxName2AlnCoreFile, $AlnCoreGenome, $c, $StableScript,
   $TaxName2AlnORFile, $TaxORF, $ORFile, $cmd, $TaxName, $d, $e, $SingleAln, $TempAln,
   $fastaAln, $sto, $ORFModel, $CoreORFstoAln, $CoreORFfastaAln, $Header, $Seq);
my(@List, @CoreGenome, @Data, @FileIntoArray, @ORFAln);
my %Alns;

$MainPath = "/home/bioinformatica/CoreGenome";
$Project = $MainPath ."/". $ProjectName;
$MainList = $Project ."/". $List;

$ORFsPath = $Project ."/". "ORFs";
$CoreGenomeFile = $Project ."/". "CoreGenome.csv";
$CoreSeqPath = $Project ."/". "CoreSeq";
$CoreAlnPath = $CoreSeqPath ."/". "AlignedCore";
$StableScript = $MainPath ."/". "src" ."/". "stable.py";

$LogFile = $Project ."/". $ProjectName . ".log";

$SeqExt = ".fasta";
$fastaAln = ".aln" . $SeqExt;
$sto = ".sto";
$HmmExt = ".hmm";

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

@List = ReadFile($MainList);
$m = scalar@List;

MakeDir($CoreAlnPath);
@CoreGenome = ReadFile($CoreGenomeFile);
$n = scalar@CoreGenome;

for ($a=0; $a<$n; $a++){
  $Row = $CoreGenome[$a];
  @Data = split(",",$Row);
  $o = scalar@Data;
  push (@FileIntoArray, [@Data]);
}

for ($b=1; $b<$n; $b++){
        $ORF = $FileIntoArray[$b]->[0];
        $CoreORF = $ORFsPath ."/". $ORF ."/". "Core-" . $ORF . $SeqExt;
        $ORFModel = $ORFsPath ."/". $ORF ."/". $ORF . $HmmExt;
        $TempAln = $ORFsPath ."/". $ORF ."/". "Core-" . $ORF . "-Temp" . $fastaAln;
        $CoreORFstoAln = $ORFsPath ."/". $ORF ."/". "Core-" . $ORF . $sto;
        $CoreORFfastaAln = $ORFsPath ."/". $ORF ."/". "Core-" . $ORF . $fastaAln;
        for ($c=1; $c<$o; $c++){
              $TaxName = $FileIntoArray[0]->[$c];
              $TaxORF = $FileIntoArray[$b]->[$c];
              
              $ORFile = $ORFsPath ."/". $ORF ."/". $TaxName ."-". $TaxORF . $SeqExt;   
              $cmd = `cat $ORFile >> "$CoreORF"`;
        }
       print "\tAligning $ORF...";
       $cmd = `muscle -in $CoreORF -out $CoreORFfastaAln -quiet`;
       #$cmd = `hmmalign $ORFModel $CoreORF > $CoreORFstoAln`;
       #$cmd = `esl-reformat fasta $CoreORFstoAln > $CoreORFfastaAln`;
       print "Done!\n";
}

for ($d=1; $d<$n; $d++){
        $ORF = $FileIntoArray[$d]->[0];
        $AlnCoreORF = $ORFsPath ."/". $ORF ."/". "Core-" . $ORF . $fastaAln;
        @ORFAln = ReadMultiFastaFile($AlnCoreORF);
        
        %Alns;
        foreach $SingleAln(@ORFAln){
                ($Header, $Seq) = ReadSeq($SingleAln);
                $Alns{$Header}= $Seq;
        }

        for ($e=1; $e<$o; $e++){
                $TaxName = $FileIntoArray[0]->[$e];
                $TaxORF = $FileIntoArray[$d]->[$e]; 
                $AlnCoreGenome = $CoreAlnPath ."/". $TaxName . $fastaAln; 
                        
                print "Adding aligned $ORF into $TaxName CoreGenome file..."; 
                        
                open (FILE, ">>$AlnCoreGenome");
                        print FILE "$Header\n$Alns{$Header}\n";
                close FILE;
                print "Done!\n";
        }
}        


        #for ($e=1; $e<$o; $e++){ 
        #        $TaxName = $FileIntoArray[0]->[$e];
        #        print "Adding aligned $ORF into $TaxName CoreGenome file...";
        #        
        #        $AlnCoreGenome = $CoreAlnPath ."/". $TaxName . $AlnExt;
        #        $SingleAln = $ORFAln[$e-1];
        #        #($Header, $Seq) = ReadSeq($SingleAln);
        #        
        #        open (FILE, ">>$AlnCoreGenome");
        #                print FILE $SingleAln;
        #        close FILE;
        #        print "Done!\n";
        #}


exit;
