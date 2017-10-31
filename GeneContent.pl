#!/usr/bin/perl -w
#################################################################################
#Scipt ContigsLength.pl                                                         #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          11 de abril de 2017                                             #
#################################################################################

use strict;
use lib '/Users/bioinformatica/Documents/CoreGenome/src/lib';
use Routines;
use List::MoreUtils qw{any};

my ($Usage, $ProjectName, $List, $CPUs, $MainPath);

$Usage = "\tUsage: CoreGenome.pl <Project Name> <List File Name> <CPUs> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List = $ARGV[1];
$CPUs = $ARGV[2];
$MainPath = $ARGV[3];

my($Project, $MainList, $ORFeomesPath, $BlastPath, $ORFsPath,
   $InitialPresenceAbsence, $PresenceAbsence, $PanGenomeSeq, $Stats, $SeqExt,
   $AlnExt, $HmmExt, $TotalPresenceAbsence, $TotalQry, $QryDb, $TotalQryIDs,
   $TotalNewORFs, $CoreGenomeSize, $Count, $Counter, $QryGenomeName, $TestingORF,
   $QryGenomeSeq, $Hmm, $ORFTemp, $cmd, $ORFpath, $Line, $ORFStoAln, $Entry, $Strand,
   $QryORFSeq, $BestHit, $Row, $Fh, $Reference, $Closest, $GeneTemp, $ReferenceORFID,
   $ClosestORFID, $CoreGenomeFile, $LogFile, $Lap, $Element, $CoreSeqsPath, $ORF,
   $Strain, $Db, $OutCore, $Gene, $ORFStoAlnFragment, $LastORFAln, $NewORFAln,
   $QryORFFastaAln, $PreviousAlnPrefix, $CurrentAlnPrefix, $CoreGenes,
   $PreviousAlnName, $FindingPreviousAln, $ID, $QryId, $NewORFId, $Summary, $NewCounter,
   $NewORF, $NewORFPath, $NewORFSeq, $NewORFHmm, $PreviousAln, $nPanGenome, $NewStrains);
my($o, $j, $i, $c, $d, $e, $p, $f, $g, $h, $q, $x);
my(@List, @PresenceAbsence, @nHMMerReport, @BestHitArray, @DataInRow,
   @LastReportColumnData, @NewReport, @PresenceAbsenceArray, @PresenceAbsenceFields,
   @SharedORFsArray, @LapArray, @CoreFile, @OrfLine, @CoreData, @Temp, @StoAln,
   @AlnData, @CoreGenome, @QryIDs, @NewORFs, @SplitPreviousAlnName, @SplitAlnPath, @AnalyzedORFs);
my $NewReport = [ ];
my $PermutationsFile = [ ];
my $Statistics = [ ];

#$MainPath = "/Users/rc/CoreGenome";
$Project = $MainPath ."/". $ProjectName;

$SeqExt = ".fasta";
$AlnExt = ".aln.fasta";
$HmmExt = ".hmm";

$MainList       = $Project ."/". $List;
$ORFeomesPath   = $MainPath ."/". "ORFeomes" ."/". "Sorted" ."/". "Filtered";
$BlastPath      = $MainPath ."/". "Blast";
$ORFsPath       = $Project."/". "ORFs";
$InitialPresenceAbsence= $Project ."/". $ProjectName . "_Initial_Presence_Absence.csv";
$PresenceAbsence= $Project ."/". $ProjectName . "_Presence_Absence.csv";
$CoreGenomeFile = $Project ."/". $ProjectName . "_CoreGenome.csv";
$LogFile        = $Project ."/". $ProjectName . ".log";
$CoreSeqsPath   = $Project ."/". "CoreSequences";
$Summary        = $Project ."/". $ProjectName . "_Summary.txt";
$PanGenomeSeq   = $Project ."/". $ProjectName . "_PanGenome" . $SeqExt;
$Stats         	= $Project ."/". $ProjectName . "_Statistics.csv";

MakeDir($CoreSeqsPath);

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

@List = ReadFile($MainList);
$TotalQry = scalar@List;
@PresenceAbsence = ReadFile($InitialPresenceAbsence);
$TotalPresenceAbsence = scalar@PresenceAbsence;
system("rm $InitialPresenceAbsence");
for ($a=0; $a<$TotalPresenceAbsence; $a++){
     $Row = $PresenceAbsence[$a];
     @PresenceAbsenceFields = split(",",$Row);
     $o = scalar@PresenceAbsenceFields;
     push (@PresenceAbsenceArray, [@PresenceAbsenceFields]);
}
for ($i=0; $i<$TotalPresenceAbsence; $i++){
        for ($j=0; $j<$o; $j++){
                if (defined($PresenceAbsenceArray[$i]->[$j])){
                        $NewReport -> [$i][$j] = $PresenceAbsenceArray[$i]->[$j];
                }else{
                        $NewReport -> [$i][$j] = "";
                }
        }
}
for ($a=1; $a<$TotalQry; $a++){
        $QryGenomeName = $List[$a];
        $QryGenomeSeq = $ORFeomesPath ."/". $QryGenomeName . ".ffn";
        $QryDb = $BlastPath ."/". $QryGenomeName;
        
        $NewReport -> [0][$a+2] = $QryGenomeName;
        
        @QryIDs = AnnotatedGenes($QryGenomeSeq); #Get all names of annotated genes from query
        $TotalQryIDs = scalar@QryIDs;
        $TotalNewORFs = $TotalQryIDs;
        
        $nPanGenome = `grep ">" -c $PanGenomeSeq`;
        
        $Counter = 0;
        @AnalyzedORFs = "";
        for ($b=1; $b<$nPanGenome+1; $b++){
                
                $Counter = sprintf "%.4d", $b;
                
                $TestingORF = "ORF" ."_". $Counter;
                $ORFpath = $ORFsPath ."/". $TestingORF;
                $Hmm = $ORFpath ."/". $TestingORF . $HmmExt;
                $ORFTemp = $ORFpath ."/". $TestingORF ."-". $QryGenomeName . ".temp";
	
		print "\n----------------Looking for $TestingORF in $QryGenomeName----------------\n";
	
		system("nhmmer -E 1e-10 --cpu $CPUs --noali --dfamtblout $ORFTemp $Hmm $QryGenomeSeq");             
                
		@nHMMerReport = ReadFile($ORFTemp);
                if(@nHMMerReport){
                        open (FILE, ">$ORFTemp");
                                print FILE $nHMMerReport[0];
                        close FILE;
                        @nHMMerReport = ReadFile($ORFTemp);
                        
                        #$CoreGenomeSize++;
			$BestHit = $nHMMerReport[0];
			$BestHit =~ s/\s^//g;
			$BestHit =~ s/\s+/,/g;
			@BestHitArray = split(",",$BestHit);
			$Entry = $BestHitArray[0];
			$Strand = $BestHitArray[8];
			$QryORFSeq = $ORFpath ."/". $QryGenomeName ."-". $Entry . $SeqExt;
                        
                        if ( any { $_ eq $Entry} @AnalyzedORFs){     
                                print "\n$Entry is already analyzed\n";    
                        }else{
                                unshift @AnalyzedORFs, $Entry;
                                
                                $NewReport -> [$b][$a+2] = $Entry;
                                
                                Extract($QryGenomeName,$QryDb,$Entry,$QryORFSeq);
                                
                                $PreviousAln = "*-". $TestingORF . $AlnExt;
                                $FindingPreviousAln = `find $ORFsPath -type f -name $PreviousAln`;
                                @SplitAlnPath = split("/",$FindingPreviousAln);
                                $PreviousAlnName = $SplitAlnPath[$#SplitAlnPath];
                                @SplitPreviousAlnName = split("-",$PreviousAlnName);
                                $PreviousAlnPrefix = $SplitPreviousAlnName[0];
                                
                                $CurrentAlnPrefix = $PreviousAlnPrefix+1;
                                $LastORFAln = $ORFpath ."/". $PreviousAlnPrefix ."-". $TestingORF . $AlnExt;
                                $NewORFAln = $ORFpath ."/". $CurrentAlnPrefix ."-". $TestingORF . $AlnExt;
                                
                                print "\tUpdating hmm profile...";
                                system("hmmalign -o $NewORFAln --trim --mapali $LastORFAln $Hmm $QryORFSeq");
                                system("hmmbuild $Hmm $NewORFAln");
                                print "Done!\n";
                                
                                #Deleting the file obtained from nhmmer 
                                print "\tCleaning...";
                                system("rm $ORFTemp $LastORFAln");
                                print "Done!\n";
                                
                                #Obataining new ORFs IDs
                                for($x=0;$x<$TotalNewORFs;$x++){
                                        $ID = $QryIDs[$x];
                                        if($ID eq $Entry){
                                                splice @QryIDs, $x, 1;
                                                $TotalNewORFs--;
                                        }
                                }
                        }
		}else{
                        print "The ORF $TestingORF is not present in $QryGenomeName\n";
                        $NewReport -> [$b][$a+2] = "";
                        
                        system("rm $ORFTemp");
                }
	}    
        print "\nProcessing New ORFs\n"; # <- Non shared genes from query file
        $TotalNewORFs = scalar@QryIDs;
        for($c=0; $c<$TotalNewORFs; $c++){
                $NewORFId = $QryIDs[$c];
                $NewCounter = sprintf "%.4d",$Counter+1+$c;
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
                for ($i=1; $i<$a+1;$i++){
                        $NewReport -> [$NewCounter][$i] = "";
                }
                $NewReport -> [$NewCounter][$a+2] = $NewORFId;
        }
        open (FILE, ">$PresenceAbsence");
        for($d=0; $d<$NewCounter+1; $d++){
                for ($e=0; $e<$a+3; $e++){
                        print FILE $NewReport -> [$d][$e], ",";
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
                       if($Count == $a+3){
                           print FILE "$Lap\n";   
                       }
                }
        close FILE;

        @CoreGenome = ReadFile($CoreGenomeFile);
        $q = scalar@CoreGenome;
        $CoreGenes = $q-1;
        
        $NewStrains = $a+1;
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
$f = scalar@CoreFile;

foreach $ORF(@CoreFile){
       @OrfLine = split(",",$ORF);
       push (@CoreData, [@OrfLine]);    
}

#Building a fasta file wirh the core genes for each genome. All genes withing core genome
#files are in the same order
print "Building CoreGenomes fasta files:\n";
for ($g=2; $g<$TotalQry+2; $g++){
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