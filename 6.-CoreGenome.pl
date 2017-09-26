#!/usr/bin/perl -w

#################################################################################
#   Programa Pre-Core-Genome                                                    #
#   Nota: Se debe ajustar la ruta de lib y de la variable $PathSeq a las que    #
#   realmente se tengan donde se instalción el programa.                        #
#                                                                               #
# Programador:   Roberto C. Torres                                              #
# Fecha:         11 de abril de 2017                                            #
#################################################################################

use strict;
use lib '/Users/Roberto/lib';
use Routines;

my ($Usage, $ProjectName, $List, $CPUs);

$Usage = "\tUsage: CoreGenome.pl <Project Name> <List File Name> <CPUs>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List = $ARGV[1];
$CPUs = $ARGV[2];

my($MainPath, $Project, $MainList, $ORFeomesPath, $BlastPath, $ORFsPath,
   $PresenceAbsence, $SharedORFs, $SeqExt, $AlnExt, $HmmExt, $n, $m, $o,
   $Count, $i, $Counter, $QryGenomeName, $TestingORF, $QryGenomeSeq, $Hmm,
   $ORFTemp, $QryGenomeDb, $cmd, $ORFpath, $Line, $ORFStoAln, $Entry, $Strand,
   $QryORFSeq, $BestHit, $Row, $Fh, $c, $d, $e, $Reference, $Closest, $GeneTemp,
   $ReferenceORFID, $ClosestORFID, $CoreGenomeFile, $LogFile, $Lap, $Element, $p,
   $CoreSeqsPath, $f, $ORF, $g, $Strain, $Db, $OutCore, $h, $Gene, $ORFStoAlnFragment,
   $LastORFAln, $NewORFAln, $QryORFFastaAln, $PreviousAlnPrefix, $CurrentAlnPrefix,
   $q, $CoreGenes, $Stat);
my(@List, @PresenceAbsence, @nHMMerReport, @BestHitArray, @DataInRow,
   @LastReportColumnData, @NewReport, @PresenceAbsenceArray, @PresenceAbsenceFields,
   @SharedORFsArray, @LapArray, @CoreFile, @OrfLine, @CoreData, @Temp, @StoAln,
   @AlnData, @CoreGenome, @QryIDs);
my $NewReport = [ ];
my $PermutationsFile = [ ];

$MainPath = "/Users/Roberto/CoreGenome";
$Project = $MainPath ."/". $ProjectName;

$MainList       = $Project ."/". $List;
$ORFeomesPath   = $MainPath ."/". "ORFeomes";
$BlastPath      = $MainPath ."/". "Blast";
$ORFsPath       = $Project."/". "ORFs";
$PresenceAbsence= $Project ."/". $ProjectName . "_Initial_Presence_Absence.csv";
$SharedORFs     = $Project ."/". $ProjectName . "_Presence_Absence.csv";
$CoreGenomeFile = $Project ."/". $ProjectName . "_CoreGenome.csv";
$LogFile        = $Project ."/". $ProjectName . ".log";
$CoreSeqsPath   = $Project ."/". "CoreSequences";
$Stat           = $Project ."/". $ProjectName . "_CoreStatistics.txt";

MakeDir($CoreSeqsPath);

$SeqExt = ".fasta";
$AlnExt = ".aln.fasta";
$HmmExt = ".hmm";

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

@List = ReadFile($MainList);
$m = scalar@List;
@PresenceAbsence = ReadFile($PresenceAbsence);
$n = scalar@PresenceAbsence;

for ($a=0; $a<$n; $a++){
     $Row = $PresenceAbsence[$a];
     @PresenceAbsenceFields = split(",",$Row);
     $o = scalar@PresenceAbsenceFields;
     push (@PresenceAbsenceArray, [@PresenceAbsenceFields]);
}

$Reference = $PresenceAbsenceArray[0]->[1];
$Closest = $PresenceAbsenceArray[0]->[2];

$NewReport -> [0][0] = "";
$NewReport -> [0][1] = $Reference;
$NewReport -> [0][2] = $Closest;

for ($b=1; $b<$n; $b++){
	$TestingORF = $PresenceAbsenceArray[$b]->[0];
	$ReferenceORFID = $PresenceAbsenceArray[$b]->[1];
	$ClosestORFID = $PresenceAbsenceArray[$b]->[2];
	
	$NewReport -> [$b][0] = $TestingORF;
	$NewReport -> [$b][1] = $ReferenceORFID;
	$NewReport -> [$b][2] = $ClosestORFID;
	
	$ORFpath = $ORFsPath ."/". $TestingORF;
	$Hmm = $ORFpath ."/". $TestingORF . $HmmExt;
	
     $Count = 1;
	for ($c=1; $c<$m; $c++){
		$QryGenomeName = $List[$c];
		$QryGenomeSeq = $ORFeomesPath ."/". $QryGenomeName . $SeqExt;
		$ORFTemp = $ORFpath ."/". $QryGenomeName ."_". $TestingORF . ".temp";
		$QryGenomeDb = $BlastPath ."/". $QryGenomeName;
	
		$NewReport -> [0][$c+2] = $QryGenomeName;
          
          @QryIDs = AnnotatedGenes($QryGenomeSeq); #Get all names of annotated genes from query
	
		print "\n----------------Looking for $TestingORF in $QryGenomeName----------------\n";
	
		system("nhmmer -E 1e-90 --cpu $CPUs --noali --dfamtblout $ORFTemp $Hmm $QryGenomeSeq");
		
		@nHMMerReport = ReadFile($ORFTemp);
		open (FILE, ">$ORFTemp");
			print FILE $nHMMerReport[0];
		close FILE;
		@nHMMerReport = ReadFile($ORFTemp);
	
		if (-z $ORFTemp){
			print "The ORF $TestingORF is not present in $QryGenomeName\n";
               $NewReport -> [$b][$c+2] = "";
			system("rm $ORFTemp");
		}else{
               $Count++;
			$BestHit = $nHMMerReport[0];
			$BestHit =~ s/\s^//g;
			$BestHit =~ s/\s+/,/g;
			@BestHitArray = split(",",$BestHit);
			$Entry = $BestHitArray[0];
			$Strand = $BestHitArray[8];
               $PreviousAlnPrefix = $Count-1;
               $CurrentAlnPrefix = $Count;
               
			$QryORFSeq = $ORFpath ."/". $QryGenomeName ."-". $Entry . $SeqExt;
               $LastORFAln = $ORFpath ."/". $PreviousAlnPrefix ."-". $TestingORF . $AlnExt;
               $NewORFAln = $ORFpath ."/". $CurrentAlnPrefix ."-". $TestingORF . $AlnExt;
	
			$NewReport -> [$b][$c+2] = $Entry;
               
               Extract($QryGenomeName,$QryGenomeDb,$Entry,$QryORFSeq);
	
			#print "\tExtracting $Entry from $QryGenomeName...";
			#system("blastdbcmd -db $QryGenomeDb -dbtype nucl -entry \"$Entry\" -out $QryORFSeq");
			#print "Done!\n";
	
			print "\tUpdating hmm profile...";
               system("hmmalign -o $NewORFAln --trim --mapali $LastORFAln $Hmm $QryORFSeq");
			system("hmmbuild $Hmm $NewORFAln");
	          print "Done!\n";
                
               #Deleting the file obtained from nhmmer 
			print "\tCleaning...";
			system("rm $ORFTemp $LastORFAln");
			print "Done!\n";
		}
	}
}

#Building a file with all genes present in each genome (...SharedORFs.csv)
open (FILE, ">>$SharedORFs");
for($d=0; $d<$n; $d++){
    for ($e=0; $e<$m+1; $e++){
        print FILE $NewReport -> [$d][$e], ",";
    }
    print FILE "\n";
}
close FILE;

#Building a file with only those genes that are shared in all genomes (...CoreGenome.csv)
@SharedORFsArray = ReadFile($SharedORFs);

open (FILE,">>$CoreGenomeFile");
print FILE "$SharedORFsArray[0]\n";
foreach $Lap(@SharedORFsArray){
       @LapArray = split(",",$Lap);
       chomp@LapArray;
       $Count = 0;
       foreach $Element(@LapArray){
              #if (defined($Element)){
              if ($Element ne ""){
                     $Count++;
              }
       }
       if($Count == $m+1){
           print FILE "$Lap\n";   
       }
}
close FILE;

@CoreGenome = ReadFile($CoreGenomeFile);
$q = scalar@CoreGenome;
$CoreGenes = $q-1;

open (FILE, ">>$Stat");
        print FILE "Project Name: $ProjectName\n";
        print FILE "Included Genomes List File: $MainList\n";
        print FILE "Number of taxa: $m\n";
        print FILE "Number of CoreGenes: $CoreGenes";
close FILE;

#Loading the core file report into an array of array
@CoreFile = ReadFile($CoreGenomeFile);
$f = scalar@CoreFile;

foreach $ORF(@CoreFile){
       @OrfLine = split(",",$ORF);
       push (@CoreData, [@OrfLine]);    
}

#Building a core fasta file for each genome. All core genes withing core genome
#files are in the same order
print "Building CoreGenomes fasta files:\n";
for ($g=1; $g<$m+1; $g++){
       $Strain = $CoreData[0]->[$g];
       $Db = $BlastPath ."/". $Strain;
       $OutCore = $CoreSeqsPath ."/". $Strain . "-CoreGenome" . $SeqExt;
       print "\t$Strain...";
       for($h=1; $h<$f; $h++){ 
              $Gene = $CoreData[$h]->[$g];
              $GeneTemp = $CoreSeqsPath ."/". $Strain ."-". $Gene . ".temp";
              $cmd = `blastdbcmd -db $Db -dbtype nucl -entry "$Gene" -out $GeneTemp`;          
              $cmd = `cat $GeneTemp >> "$OutCore"`;   
              $cmd = `rm $GeneTemp`;
       }
       print "Done!\n";
}
exit;

#################################################################################
sub AnnotatedGenes{
        my ($File) = @_;
        my $cmd = `grep ">" $File`;
           $cmd =~ s/>//g;
           $cmd =~ s/\h//g;
        my @Data = split('\n',$cmd);
        return @Data;
}

sub Extract{
        my ($Qry, $DataBase,$Entry,$OutSeq, $null) = @_;
        print "\tExtracting ORF from $Qry...";	
        $cmd = `blastdbcmd -db $DataBase -dbtype nucl -entry "$Entry" -out $OutSeq`;
        print "Done!\n";
}