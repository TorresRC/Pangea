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
use FindBin;
use lib "$FindBin::Bin/lib";
use Routines;
use List::MoreUtils qw{any};

my ($Usage, $ProjectName, $List, $CPUs, $MainPath, $Recovery, $Add, $AddList);

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

my($Project, $SeqExt, $AlnExt, $stoExt, $HmmExt, $SharedORFs, $CoreGenomeFile,
   $TempFile, $MainList, $ORFeomesPath, $BlastPath, $ORFsPath, $PresenceAbsence,
   $LogFile, $CoreSeqsPath, $TotalQry, $Qry, $QryIn, $QryDb,
   $LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, $Line, $Count, $Element,
   $LinesOnCoreGenome, $ColumnsOnCoreGenome, $KnownORF, $ORFpath, $Hmm, $ORFTemp,
   $BestHit, $Entry, $Strand, $QryORFSeq, $PreviousAln, $CheckAln,
   $FindingPreviousAln, $PreviousAlnName, $PreviousAlnPrefix, $CurrentAlnPrefix,
   $LastORFAln, $NewORFAln, $Strain, $Db, $OutCore, $Gene, $cmd);
my($i,$j,$k);
my(@List, @PresenceAbsenceArray, @PresenceAbsenceMatrix, @Fields,
   @CoreGenomeMatrix, @AnalyzedORFs, @nHMMerReport, @BestHitArray, @SplitAlnPath,
   @SplitPreviousAlnName);
my $NewReport = [ ];

$Project = $MainPath ."/". $ProjectName;

$SeqExt = ".fasta";
$AlnExt = ".aln" . $SeqExt;
$stoExt = ".sto";
$HmmExt = ".hmm";

#################################################################################

$MainList        = $Project ."/". $List;
$ORFeomesPath    = $Project ."/". "ORFeomes" ."/". "Sorted" ."/". "Filtered";
$BlastPath       = $Project ."/". "Blast";
$ORFsPath        = $Project."/". "ORFs";
$CoreSeqsPath    = $Project ."/". "CoreSequences";
$PresenceAbsence = $Project ."/". $ProjectName . "_Initial_Presence_Absence.csv";
$CoreGenomeFile  = $Project ."/". $ProjectName . "_CoreGenome.csv";
$TempFile        = $Project ."/". "temp";

$LogFile         = $Project ."/". $ProjectName . ".log";

#################################################################################

open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

@List = ReadFile($MainList);
$TotalQry = scalar@List;

for ($i=2; $i<$TotalQry; $i++){
        $Qry = $List[$i];
        $QryIn = $ORFeomesPath ."/". $Qry . ".ffn";
        $QryDb = $BlastPath ."/". $Qry;

        @PresenceAbsenceArray = ReadFile($PresenceAbsence);
        ($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, @PresenceAbsenceMatrix) = Matrix($PresenceAbsence);

        open (FILE,">$PresenceAbsence");
                #print FILE "$PresenceAbsenceArray[0]\n";
        foreach $Line(@PresenceAbsenceArray){
                @Fields = split(",",$Line);
                chomp@Fields;
                $Count = 0;
                foreach $Element(@Fields){
                        if ($Element ne ""){
                                $Count++;
                        }
                }
                if($Count == $ColumnsOnPresenceAbsence){
                        print FILE "$Line\n";
                }
        }
        close FILE;

        ($LinesOnCoreGenome, $ColumnsOnCoreGenome, @CoreGenomeMatrix) = Matrix($PresenceAbsence);
        
        #Initializing the new presence absence report with the already analized data
        for ($j=0; $j<$LinesOnCoreGenome; $j++){
                for ($k=0; $k<$ColumnsOnCoreGenome; $k++){
                        if (defined($CoreGenomeMatrix[$j]->[$k])){
                                $NewReport -> [$j][$k] = $CoreGenomeMatrix[$j]->[$k];
                        }else{
                                $NewReport -> [$j][$k] = "";
                        }
                }
        }
        
        $NewReport -> [0][$i+1] = $Qry;
        
        @AnalyzedORFs = "";
        for ($j=1; $j<$LinesOnCoreGenome; $j++){
                $KnownORF = $CoreGenomeMatrix[$j]->[0];
                $ORFpath = $ORFsPath ."/". $KnownORF;
                $Hmm = $ORFpath ."/". $KnownORF . $HmmExt;
                $ORFTemp = $ORFpath ."/". $KnownORF ."-". $Qry . ".temp";
                
                print "\n----------------Looking for $KnownORF in $Qry----------------\n";
                
                system("nhmmer -E 1e-10 --cpu $CPUs --noali --dfamtblout $ORFTemp $Hmm $QryIn");
                
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
                        $QryORFSeq = $ORFpath ."/". $Qry ."-". $Entry . $SeqExt;
                        
                        #If the Entry is not already analyzed, extract it and incloude it into the hmm 
                        if ( any { $_ eq $Entry} @AnalyzedORFs){  
                                print "\n$Entry is already analyzed\n";    
                        }else{
                                unshift @AnalyzedORFs, $Entry;
                                
                                $NewReport -> [$j][$i+1] = $Entry;
                                
                                if (not -e $QryORFSeq){
                                        Extract($Qry,$QryDb,$Entry,$QryORFSeq);
                                }
                                        
                                $PreviousAln = $ORFpath ."/*-". $KnownORF . $AlnExt;
                                
                                # Checking if this entry is already included in the current ORF alignment
                                $CheckAln = `grep $Entry $PreviousAln`;
                                if ($CheckAln eq ""){
                                        $FindingPreviousAln = `find $PreviousAln`;
                                        
                                        @SplitAlnPath = split("/",$FindingPreviousAln);
                                        $PreviousAlnName = $SplitAlnPath[$#SplitAlnPath];
                                        @SplitPreviousAlnName = split("-",$PreviousAlnName);
                                        $PreviousAlnPrefix = $SplitPreviousAlnName[0];
                                        
                                        $CurrentAlnPrefix = $PreviousAlnPrefix+1;
                                        $LastORFAln = $ORFpath ."/". $PreviousAlnPrefix ."-". $KnownORF . $AlnExt;
                                        $NewORFAln = $ORFpath ."/". $CurrentAlnPrefix ."-". $KnownORF . $AlnExt;
                                        
                                        print "\tUpdating hmm previous align...";
                                        system("hmmalign -o $NewORFAln --trim --mapali $LastORFAln $Hmm $QryORFSeq");
                                        system("hmmbuild $Hmm $NewORFAln");
                                        print "Done!\n";
                                        
                                        system("rm $LastORFAln");
                                }
                        }
                }else{
                        print "The ORF $KnownORF is not present in $Qry\n";
                        $NewReport -> [$j][$i+1] = "";
                }
                system("rm $ORFTemp");
        }
        open (FILE, ">$PresenceAbsence");
                for($j=0; $j<$LinesOnCoreGenome; $j++){
                        for ($k=0; $k<$i+2; $k++){
                                if (defined ($NewReport -> [$j][$k])){
                                        print FILE $NewReport -> [$j][$k], ",";      
                                }else{
                                        $NewReport -> [$j][$k] = "";
                                        print FILE $NewReport -> [$j][$k], 
                                }
                        }
                        print FILE "\n";
                }
        close FILE;
}

open (FILE,">$PresenceAbsence");
        foreach $Line(@PresenceAbsenceArray){
                @Fields = split(",",$Line);
                chomp@Fields;
                $Count = 0;
                foreach $Element(@Fields){
                        if ($Element ne ""){
                                $Count++;
                        }
                }
                if($Count == $ColumnsOnPresenceAbsence){
                        print FILE "$Line\n";
                }
        }
close FILE;

system("mv $PresenceAbsence $CoreGenomeFile");
($LinesOnCoreGenome, $ColumnsOnCoreGenome, @CoreGenomeMatrix) = Matrix($CoreGenomeFile);
MakeDir($CoreSeqsPath);

for ($i=1; $i<$ColumnsOnCoreGenome; $i++){
       $Strain = $CoreGenomeMatrix[0]->[$i];
       $Db = $BlastPath ."/". $Strain;
       $OutCore = $CoreSeqsPath ."/". $Strain . "-CoreGenome" . $SeqExt;
       print "Building $Strain CoreGenome file...";
       for($j=1; $j<$LinesOnCoreGenome; $j++){ 
              $Gene = $CoreGenomeMatrix[$j]->[$i];
              $cmd = `blastdbcmd -db $Db -dbtype nucl -entry "$Gene" -out $TempFile`;          
              $cmd = `cat $TempFile >> "$OutCore"`;   
              $cmd = `rm $TempFile`;
       }
       print "Done!\n";
}
exit;