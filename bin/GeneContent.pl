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
use List::MoreUtils qw{any};
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $List, $CPUs, $MainPath, $Recovery, $Add, $AddList, $eValue);

$Usage = "\tUsage: CoreGenome.pl <Project Name> <List File Name> <CPUs> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List        = $ARGV[1];
$CPUs        = $ARGV[2];
$MainPath    = $ARGV[3];
$Recovery    = $ARGV[4];
$Add         = $ARGV[5];
$AddList     = $ARGV[6];
$eValue      = $ARGV[7];

my($Project, $MainList, $ORFeomesPath, $BlastPath, $ORFsPath,
   $InitialPresenceAbsence, $PresenceAbsence, $PanGenomeSeq, $Stats, $SeqExt,
   $AlnExt, $HmmExt, $TotalQry, $QryDb, $TotalQryIDs, $TotalNewORFs,
   $CoreGenomeSize, $Count, $Counter, $QryGenomeName, $TestingORF, $QryGenomeSeq,
   $Hmm, $ORFTemp, $cmd, $ORFpath, $Line, $ORFStoAln, $Entry, $Strand,
   $QryORFSeq, $BestHit, $Row, $Fh, $Reference, $Closest, $GeneTemp,
   $ReferenceORFID, $ClosestORFID, $CoreGenomeFile, $LogFile, $Lap, $Element,
   $CoreSeqsPath, $ORF, $Strain, $Db, $OutCore, $Gene, $ORFStoAlnFragment,
   $LastORFAln, $NewORFAln, $QryORFFastaAln, $PreviousAlnPrefix,
   $CurrentAlnPrefix, $CoreGenes, $PreviousAlnName, $FindingPreviousAln, $ID,
   $QryId, $NewORFId, $Summary, $NewCounter, $NewORF, $NewORFPath, $NewORFSeq,
   $NewORFHmm, $PreviousAln, $nPanGenome, $NewStrains, $Progress, $CheckAln,
   $DbHeaders);
my($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, $LinesOnCoreGenome,
   $nCore, $i, $j, $k);
my(@List, @PresenceAbsence, @nHMMerReport, @BestHitArray, @DataInRow,
   @LastReportColumnData, @NewReport, @PresenceAbsenceArray,
   @PresenceAbsenceFields, @SharedORFsArray, @LapArray, @CoreFile, @OrfLine,
   @CoreData, @Temp, @StoAln, @AlnData, @CoreGenome, @QryIDs, @NewORFs,
   @SplitPreviousAlnName, @SplitAlnPath, @AnalyzedORFs, @Progress);
my $NewReport = [ ];
my $PermutationsFile = [ ];
my $Statistics = [ ];

$Project = $MainPath ."/". $ProjectName;

$SeqExt = ".fasta";
$AlnExt = ".aln.fasta";
$HmmExt = ".hmm";

$MainList               = $Project ."/". $List;
$ORFeomesPath           = $Project ."/". "ORFeomes" ."/". "Sorted" ."/". "Filtered";
$BlastPath              = $Project ."/". "Blast";
$ORFsPath               = $Project."/". "ORFs";
$InitialPresenceAbsence = $Project ."/". $ProjectName . "_Initial_Presence_Absence.csv";
$PresenceAbsence        = $Project ."/". $ProjectName . "_Presence_Absence.csv";
$CoreGenomeFile         = $Project ."/". $ProjectName . "_CoreGenome.csv";
$CoreSeqsPath           = $Project ."/". "CoreSequences";
$Summary                = $Project ."/". $ProjectName . "_Summary.txt";
$PanGenomeSeq           = $Project ."/". $ProjectName . "_PanGenome" . $SeqExt;
$Stats                  = $Project ."/". $ProjectName . "_Statistics.csv";

$LogFile                = $Project ."/". $ProjectName . ".log";

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

#Starting from scratch or recovering evaluation 
if($Recovery == "0"){
        $Progress = 1;
        @PresenceAbsence = ReadFile($InitialPresenceAbsence);
        #system("rm $InitialPresenceAbsence");
}elsif($Recovery == "1"){
        @Progress = ReadFile($Stats);
        $Progress = scalar@Progress-2;
        @PresenceAbsence = ReadFile($PresenceAbsence);
}
$LinesOnPresenceAbsence = scalar@PresenceAbsence;

#Loading query names
#if($Add == "0"){
        @List = ReadFile($MainList); 
#}elsif($Add == "1"){
#        @List = ReadFile($AddList);
#}
$TotalQry = scalar@List;

#Loading intitial comparison or previous job into an array of arrays
for ($i=0; $i<$LinesOnPresenceAbsence; $i++){
        $Line = $PresenceAbsence[$i];
        @PresenceAbsenceFields = split(",",$Line);
        $ColumnsOnPresenceAbsence = scalar@PresenceAbsenceFields;
        push (@PresenceAbsenceArray, [@PresenceAbsenceFields]);
}

#Initializing the new presence absence report with the already analized data 
for ($i=0; $i<$LinesOnPresenceAbsence; $i++){
        for ($j=0; $j<$ColumnsOnPresenceAbsence; $j++){
                if (defined($PresenceAbsenceArray[$i]->[$j])){
                        $NewReport -> [$i][$j] = $PresenceAbsenceArray[$i]->[$j];
                }else{
                        $NewReport -> [$i][$j] = "";
                }
        }
}
        

#Job for each query
for ($i=$Progress; $i<$TotalQry; $i++){
        $QryGenomeName = $List[$i]; 
        $QryGenomeSeq = $ORFeomesPath ."/". $QryGenomeName . ".ffn";
        $QryDb = $BlastPath ."/". $QryGenomeName;
        $DbHeaders = $QryDb . ".nhr";
        
        if (not -e $DbHeaders){
                $cmd = `makeblastdb -in $QryGenomeSeq -dbtype nucl -parse_seqids -out $QryDb`;
        }
        
        $NewReport -> [0][$i+2] = $QryGenomeName;
        
        @QryIDs = AnnotatedGenes($QryGenomeSeq); # Get all names of query annotated genes
        $TotalQryIDs = scalar@QryIDs;
        $TotalNewORFs = $TotalQryIDs; # Starting with all query genes
        
        $nPanGenome = `grep ">" -c $PanGenomeSeq`; # Get the number of genes on the PanGenome
        
        #ORF level
        # Processing all shared genes
        $Counter = 0;
        @AnalyzedORFs = "";
        for ($j=1; $j<$nPanGenome+1; $j++){
                
                $Counter = sprintf "%.4d", $j; #Using the numeber of PanGenome's genes as index
                $TestingORF = "ORF" ."_". $Counter;
                $ORFpath = $ORFsPath ."/". $TestingORF;
                $Hmm = $ORFpath ."/". $TestingORF . $HmmExt;
                $ORFTemp = $ORFpath ."/". $TestingORF ."-". $QryGenomeName . ".temp";
                
                print "\n----------------Looking for $TestingORF in $QryGenomeName----------------\n";
                
                system("nhmmer -E $eValue --cpu $CPUs --noali --dfamtblout $ORFTemp $Hmm $QryGenomeSeq");
                
                #Analyzing the hmm results
                @nHMMerReport = ReadFile($ORFTemp);
                if(@nHMMerReport){
                        open (FILE, ">$ORFTemp");
                                print FILE $nHMMerReport[0];
                        close FILE;
                        @nHMMerReport = ReadFile($ORFTemp);
                        
                        #Entry level
                        #$CoreGenomeSize++;
                        $BestHit = $nHMMerReport[0];
                        $BestHit =~ s/\s^//g;
                        $BestHit =~ s/\s+/,/g;
                        @BestHitArray = split(",",$BestHit);
                        $Entry = $BestHitArray[0];
                        $Strand = $BestHitArray[8];
                        $QryORFSeq = $ORFpath ."/". $QryGenomeName ."-". $Entry . $SeqExt;
                        
                        #If the Entry is not already analyzed, extract it and incloude it into the hmm 
                        if ( any { $_ eq $Entry} @AnalyzedORFs){  
                                print "\n$Entry is already analyzed\n";    
                        }else{
                                unshift @AnalyzedORFs, $Entry;
                                
                                $NewReport -> [$j][$i+2] = $Entry;
                                
                                if (not -e $QryORFSeq){
                                        Extract($QryGenomeName,$QryDb,$Entry,$QryORFSeq);
                                }
                                        
                                $PreviousAln = $ORFpath ."/*-". $TestingORF . $AlnExt;
                                
                                # Checking if this entry is already included in the current ORF alignment
                                $CheckAln = `grep $Entry $PreviousAln`;
                                if ($CheckAln eq ""){
                                        $FindingPreviousAln = `find $PreviousAln`;
                                        
                                        @SplitAlnPath = split("/",$FindingPreviousAln);
                                        $PreviousAlnName = $SplitAlnPath[$#SplitAlnPath];
                                        @SplitPreviousAlnName = split("-",$PreviousAlnName);
                                        $PreviousAlnPrefix = $SplitPreviousAlnName[0];
                                        
                                        $CurrentAlnPrefix = $PreviousAlnPrefix+1;
                                        $LastORFAln = $ORFpath ."/". $PreviousAlnPrefix ."-". $TestingORF . $AlnExt;
                                        $NewORFAln = $ORFpath ."/". $CurrentAlnPrefix ."-". $TestingORF . $AlnExt;
                                        
                                        print "\tUpdating hmm previous alignment...";
                                        system("hmmalign -o $NewORFAln --trim --mapali $LastORFAln $Hmm $QryORFSeq");
                                        system("hmmbuild $Hmm $NewORFAln");
                                        print "Done!\n";
                                        
                                        system("rm $LastORFAln");
                                }
                                
                                # Delete the Entry of the TotalNewORFs colection and keeping only the non shared genes
                                for($k=0;$k<$TotalNewORFs;$k++){
                                        $ID = $QryIDs[$k];
                                        if($ID eq $Entry){
                                                splice @QryIDs, $k, 1;
                                                $TotalNewORFs--;
                                        }
                                }
                        }
                }else{
                        print "The ORF $TestingORF is not present in $QryGenomeName\n";
                        $NewReport -> [$j][$i+2] = "";
                }
                system("rm $ORFTemp");
        }    
        
        
        # Processing all non shared genes from query file
        print "\nProcessing New ORFs\n";
        $TotalNewORFs = scalar@QryIDs; # Number of remaining genes 
        for($j=0; $j<$TotalNewORFs; $j++){
                $NewORFId = $QryIDs[$j];
                $NewCounter = sprintf "%.4d",$Counter+1+$j;
                $NewORF = "ORF" ."_". $NewCounter;
                $NewORFPath = $ORFsPath ."/". $NewORF;
                $NewORFSeq = $NewORFPath ."/". $QryGenomeName ."-". $NewORFId . $SeqExt;
                $NewORFHmm = $NewORFPath ."/". $NewORF . $HmmExt;
                $NewORFAln = $NewORFPath ."/". "1-" . $NewORF . $AlnExt;
                
                print "\nProcessing ORF $NewCounter: \n";
                MakeDir($NewORFPath);
                Extract($QryGenomeName,$QryDb,$NewORFId,$NewORFSeq);
                $cmd = `cp $NewORFSeq $NewORFAln`;
                HMM($CPUs,$NewORFHmm,$NewORFAln);
                        
                print "\tAdding ORF $NewCounter to PanGenome...";
                $cmd = `blastdbcmd -db $QryDb -dbtype nucl -entry "$NewORFId" >> $PanGenomeSeq`;
                print "Done!\n";
                        
                $NewReport -> [$NewCounter][0] = $NewORF;
                
                for ($k=1; $k<$i+1;$k++){
                        $NewReport -> [$NewCounter][$k] = "";
                }
                $NewReport -> [$NewCounter][$i+2] = $NewORFId;
        }

        # Building a Presence/Absence file
        open (FILE, ">$PresenceAbsence");
        for($j=0; $j<$NewCounter+1; $j++){
                for ($k=0; $k<$i+3; $k++){
                        if (defined ($NewReport -> [$j][$k])){
                                print FILE $NewReport -> [$j][$k], ",";      
                        }else{
                                $NewReport -> [$j][$k] = "";
                                print FILE $NewReport -> [$j][$k], ",";
                        }
                }
                print FILE "\n";
        }
        close FILE;
        
        #Building a file with only those genes that are shared in all genomes (Core-Genome table)
        @SharedORFsArray = ReadFile($PresenceAbsence);  
        open (FILE,">$CoreGenomeFile");
                #print FILE "$SharedORFsArray[0]\n";
                foreach $Lap (@SharedORFsArray){
                        @LapArray = split(",",$Lap);
                        chomp@LapArray;
                        $Count = 0;
                        foreach $Element(@LapArray){
                                if ($Element ne ""){
                                        $Count++;
                                }
                        }
                        if($Count == $i+3){
                                print FILE "$Lap\n";   
                        }
                }
        close FILE;
        
        @CoreGenome = ReadFile($CoreGenomeFile);
        $LinesOnCoreGenome = scalar@CoreGenome;
        $CoreGenes = $LinesOnCoreGenome-1;
                
        $NewStrains = $i+1;
        open (FILE, ">>$Stats");
        #  "Number Of New Strains,Analyzed Strain,Core Genome,Pan Genome,New Genes\n";
        print FILE "$NewStrains,$QryGenomeName,$CoreGenes,$NewCounter,$TotalNewORFs\n";
        close FILE;
}


open (FILE, ">>$Summary");
        print FILE "Core-Genome Size: $CoreGenes\nPan-Genome Size: $NewCounter";
close FILE;


#Loading the core file report into an array of array
@CoreFile = ReadFile($CoreGenomeFile);
$nCore = scalar@CoreFile;

MakeDir($CoreSeqsPath);

foreach $ORF(@CoreFile){
       @OrfLine = split(",",$ORF);
       push (@CoreData, [@OrfLine]);    
}
#Building a fasta file wirh the core genes for each genome. All genes withing core genome
#files are in the same order
print "Building CoreGenomes fasta files:\n";
for ($i=2; $i<$TotalQry+2; $i++){
       $Strain = $CoreData[0]->[$i];
       $Db = $BlastPath ."/". $Strain;
       $OutCore = $CoreSeqsPath ."/". $Strain . "-CoreGenome" . $SeqExt;
       print "\t$Strain...";
       for($j=1; $j<$nCore; $j++){ 
              $Gene = $CoreData[$j]->[$i];
              $GeneTemp = $CoreSeqsPath ."/". $Strain ."-". $Gene . ".temp";
              $cmd = `blastdbcmd -db $Db -dbtype nucl -entry "$Gene" -out $GeneTemp`;          
              $cmd = `cat $GeneTemp >> "$OutCore"`;   
              $cmd = `rm $GeneTemp`;
       }
       print "Done!\n";
}
exit;