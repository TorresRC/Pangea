#!/usr/bin/perl -w

#################################################################################
#   Programa CoreGenome                                                         #
#                                                                               #
# Programador:   Roberto C. Torres                                              #
# Fecha:         11 de abril de 2017                                            #
#################################################################################
use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;
use List::Util qw(any max);
use List::MoreUtils qw(uniq first_index);

my ($Usage, $ProjectName, $List, $CPUs, $Recovery, $Add, $AddList,
    $OutPath);

$Usage = "\tUsage: CoreGenome.pl <Project Name> <List File Name> <CPUs> <Out_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List        = $ARGV[1];
$CPUs        = $ARGV[2];
$OutPath     = $ARGV[3];

my($SeqExt, $AlnExt, $stoExt, $HmmExt, $SharedORFs, $CoreGenomeFile,
   $TempFile, $ORFeomesPath, $BlastPath, $ORFsPath, $PresenceAbsence,
   $LogFile, $CoreSeqsPath, $TotalQry, $Qry, $QryIn, $QryDb,
   $LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, $Line, $Count, $Element,
   $LinesOnCoreGenome, $ColumnsOnCoreGenome, $KnownORF, $Hmm, $ORFTemp,
   $BestHit, $Entry, $Strand, $QryORFSeq, $PreviousAln, $CheckAln,
   $FindingPreviousAln, $PreviousAlnName, $PreviousAlnPrefix, $CurrentAlnPrefix,
   $LastORFAln, $NewORFAln, $Strain, $Db, $OutCore, $Gene, $cmd, $CoreGenomeHmmDb,
   $ScanTempFile, $ORF, $SharedORF, $MaxScore, $Score, $QryEntry, $ORFpath,
   $ORFhmm, $EntrySeq);
my($i,$j,$k, $nCoreGenes);
my(@List, @PresenceAbsenceArray, @PresenceAbsenceMatrix, @Fields, @CoreGenes,
   @CoreGenomeMatrix, @AnalyzedORFs, @nHMMerReport, @BestHitArray, @SplitAlnPath,
   @SplitPreviousAlnName, @SharedORFs, @HmmScanFile, @HitData, @TargetName,
   @ORFName);
my(%Qry, %Entry, %Scores);
my $NewReport = [ ];

$SeqExt = ".fasta";
$AlnExt = ".aln" . $SeqExt;
$stoExt = ".sto";
$HmmExt = ".hmm";

$ORFeomesPath    = $OutPath ."/". "ORFeomes";
$BlastPath       = $OutPath ."/". "Blast";
$ORFsPath        = $OutPath ."/". "ORFs";
$CoreSeqsPath    = $OutPath ."/". "CoreSequences";
$PresenceAbsence = $OutPath ."/". $ProjectName . "_Initial_Presence_Absence.csv";
$CoreGenomeFile  = $OutPath ."/". $ProjectName . "_CoreGenome.csv";
$TempFile        = $OutPath ."/". "temp";
$CoreGenomeHmmDb = $CoreSeqsPath ."/". "CoreGenomeDb.hmm";

$LogFile         = $OutPath ."/". $ProjectName . ".log";
open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

MakeDir($CoreSeqsPath);

@List = ReadFile($List);
$TotalQry = scalar@List;

#@PresenceAbsenceArray = ReadFile($PresenceAbsence);
($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, @PresenceAbsenceMatrix) = Matrix($PresenceAbsence);

#Initializing the new presence absence report with the already analized data
for ($i=0; $i<$LinesOnPresenceAbsence; $i++){
   $ORF = $PresenceAbsenceMatrix[$i]->[0];
   #$NewReport -> [$i][0] = $PresenceAbsenceMatrix[$i]->[0];
   $Count = 0;
   for ($j=0; $j<$ColumnsOnPresenceAbsence; $j++){
      $Qry = $PresenceAbsenceMatrix[0]->[$j];
      #$NewReport -> [0][$j] = $PresenceAbsenceMatrix[0]->[$j];
      if (defined($PresenceAbsenceMatrix[$i]->[$j])){
         #$NewReport -> [$i][$j] = $PresenceAbsenceMatrix[$i]->[$j];
         $Count++;
         if ($i > 0){
            if ($j > 0){
               $Entry{$ORF}{$Qry} = $PresenceAbsenceMatrix[$i]->[$j];
            }
         }         
      }
   }
   if ($j > 0){
      $Gene = $PresenceAbsenceMatrix[$i][0];
      if($Count == $ColumnsOnPresenceAbsence){
         push @CoreGenes, $Gene;
      }
   }
}

foreach $Gene (@CoreGenes){
   $Hmm = $ORFsPath ."/". $Gene ."/". $Gene . $HmmExt;
   system("cat $Hmm >> $CoreGenomeHmmDb");
}
$nCoreGenes = scalar@CoreGenes;
system ("hmmpress $CoreGenomeHmmDb");

for ($i=1; $i<$TotalQry; $i++){
   $Qry   = $List[$i];
   $QryIn = $ORFeomesPath ."/". $Qry . ".ffn";
   $QryDb = $BlastPath ."/". $Qry;
   #$NewReport -> [0][$i+1] = $Qry;
   $ScanTempFile = $CoreSeqsPath ."/". "$Qry" . "_ScanTemp";
   
   system ("hmmscan -E 1e-3 --cpu $CPUs --noali --tblout $ScanTempFile $CoreGenomeHmmDb $QryIn");
   
   @SharedORFs="";
   @HmmScanFile = ReadFile($ScanTempFile);
   if (@HmmScanFile){
      foreach $Line(@HmmScanFile){
         $Line =~ s/\s^//g;
         $Line =~ s/\s+/,/g;
         @HitData = split(",",$Line);
         @TargetName = split("-",$HitData[0]);
         @ORFName = split('\.',$TargetName[1]);
         $ORF = $ORFName[0];
         $QryEntry = $HitData[2];
         $Score = $HitData[5];
         #if (any {$_ eq $ORF} @SharedORFs){
         #}else{
           push @SharedORFs, $ORF;
         #}
         $Qry{$ORF}{$Score}{$Qry} = $QryEntry;
         push (@{$Scores{$ORF}{$Qry}}, $Score);   
      }
   }
   shift@SharedORFs;
   @SharedORFs = uniq(@SharedORFs);
      
   system("rm $ScanTempFile $CoreGenomeHmmDb*");
   foreach $ORF (@SharedORFs){
      $MaxScore = max(@{$Scores{$ORF}{$Qry}});
      $Entry = $Qry{$ORF}{$MaxScore}{$Qry};
      $Entry{$ORF}{$Qry} = $Entry;
      $Hmm = $ORFsPath ."/". $ORF ."/". $ORF . $HmmExt;
      system("cat $Hmm >> $CoreGenomeHmmDb");
   }
   system ("hmmpress $CoreGenomeHmmDb");
}

$NewReport -> [0][0] = "ORF";

@SharedORFs = sort@SharedORFs;
$nCoreGenes = scalar@SharedORFs;

for ($i=0; $i<$nCoreGenes; $i++){
   $ORF = $SharedORFs[$i];
   $NewReport -> [$i+1][0] = $ORF;
   for ($j=0; $j<$TotalQry; $j++){
      $Qry = $List[$j];
      $NewReport -> [0][$j+1] = $Qry;
   }
}

for ($i=0; $i<$nCoreGenes; $i++){
   for ($j=0; $j<$TotalQry+1; $j++){
      $NewReport -> [$i+1][$j+1] = $Entry{$SharedORFs[$i]}{$List[$j]};
   }
}

open (FILE, ">$CoreGenomeFile");
for($i=0; $i<$nCoreGenes; $i++){
   for ($j=0; $j<$TotalQry+1; $j++){
      if (defined ($NewReport -> [$i][$j])){
         print FILE $NewReport -> [$i][$j], ",";
      }else{
         $NewReport -> [$i][$j] = "";
         print FILE $NewReport -> [$i][$j], ",";
      }
   }
        print FILE "\n";
}
close FILE;

for ($i=0; $i<$TotalQry; $i++){
       $Strain = $List[$i];
       $Db = $BlastPath ."/". $Strain;
       $OutCore = $CoreSeqsPath ."/". $Strain . "-CoreGenome" . $SeqExt;
       print "Building $Strain CoreGenome file...";
       for($j=0; $j<scalar@SharedORFs; $j++){
         $ORF = $SharedORFs[$j];
         $Entry = $Entry{$ORF}{$Strain};
         
         $ORFpath = $ORFsPath ."/". $ORF;
         $LastORFAln = $ORFpath ."/". $i+1 ."_". $ORF . ".aln.fasta";
         $NewORFAln = $ORFpath ."/". $i+2 ."_". $ORF .".aln.fasta";
         $ORFhmm = $ORFpath ."/". $ORF . ".hmm";
         $EntrySeq = $ORFpath ."/". $Entry . ".fasta";
         
         $cmd = `blastdbcmd -db $Db -dbtype nucl -entry "$Entry" -out $EntrySeq`;
         $cmd = `cat $EntrySeq >> "$OutCore"`;
         
         system("hmmalign -o $NewORFAln --mapali $LastORFAln $ORFhmm $EntrySeq");
         system("hmmbuild $ORFhmm $NewORFAln");
         system("rm $LastORFAln");
         
       }
       print "Done!\n";
}

exit;

#open (FILE, ">$PresenceAbsence");
#        for($j=0; $j<$LinesOnCoreGenome; $j++){
#                for ($k=0; $k<$i+2; $k++){
#                        if (defined ($NewReport -> [$j][$k])){
#                                print FILE $NewReport -> [$j][$k], ",";      
#                        }else{
#                                $NewReport -> [$j][$k] = "";
#                                print FILE $NewReport -> [$j][$k], 
#                        }
#                }
#                print FILE "\n";
#        }
#close FILE;
 
        
        @AnalyzedORFs = "";
        #for ($j=0; $j<$LinesOnCoreGenome; $j++){
        #        $KnownORF = $CoreGenomeMatrix[$j]->[0];
        #        $ORFpath = $ORFsPath ."/". $KnownORF;
        #        $Hmm = $ORFpath ."/". $KnownORF . $HmmExt;
        #        $ORFTemp = $ORFpath ."/". $KnownORF ."-". $Qry . ".temp";
        #        
        #        print "\n----------------Looking for $KnownORF in $Qry----------------\n";
        #        
        #        system("nhmmer -E 1e-10 --cpu $CPUs --noali --dfamtblout $ORFTemp $Hmm $QryIn");
        #        
        #        #Analyzing the hmm results
        #        @nHMMerReport = ReadFile($ORFTemp);
        #        if(@nHMMerReport){
        #                open (FILE, ">$ORFTemp");
        #                        print FILE $nHMMerReport[0];
        #                close FILE;
        #                @nHMMerReport = ReadFile($ORFTemp);
        #                
        #                #Entry level
        #                #$CoreGenomeSize++;
        #                $BestHit = $nHMMerReport[0];
        #                $BestHit =~ s/\s^//g;
        #                $BestHit =~ s/\s+/,/g;
        #                @BestHitArray = split(",",$BestHit);
        #                $Entry = $BestHitArray[0];
        #                $Strand = $BestHitArray[8];
        #                $QryORFSeq = $ORFpath ."/". $Qry ."-". $Entry . $SeqExt;
        #                
        #                #If the Entry is not already analyzed, extract it and incloude it into the hmm 
        #                if ( any { $_ eq $Entry} @AnalyzedORFs){  
        #                        print "\n$Entry is already analyzed\n";    
        #                }else{
        #                        unshift @AnalyzedORFs, $Entry;
        #                        
        #                        $NewReport -> [$j][$i+1] = $Entry;
        #                        
        #                        if (not -e $QryORFSeq){
        #                                Extract($Qry,$QryDb,$Entry,$QryORFSeq);
        #                        }
        #                                
        #                        $PreviousAln = $ORFpath ."/*-". $KnownORF . $AlnExt;
        #                        
        #                        # Checking if this entry is already included in the current ORF alignment
        #                        $CheckAln = `grep $Entry $PreviousAln`;
        #                        if ($CheckAln eq ""){
        #                                $FindingPreviousAln = `find $PreviousAln`;
        #                                
        #                                @SplitAlnPath = split("/",$FindingPreviousAln);
        #                                $PreviousAlnName = $SplitAlnPath[$#SplitAlnPath];
        #                                @SplitPreviousAlnName = split("-",$PreviousAlnName);
        #                                $PreviousAlnPrefix = $SplitPreviousAlnName[0];
        #                                
        #                                $CurrentAlnPrefix = $PreviousAlnPrefix+1;
        #                                $LastORFAln = $ORFpath ."/". $PreviousAlnPrefix ."-". $KnownORF . $AlnExt;
        #                                $NewORFAln = $ORFpath ."/". $CurrentAlnPrefix ."-". $KnownORF . $AlnExt;
        #                                
        #                                print "\tUpdating hmm previous align...";
        #                                system("hmmalign -o $NewORFAln --trim --mapali $LastORFAln $Hmm $QryORFSeq");
        #                                system("hmmbuild $Hmm $NewORFAln");
        #                                print "Done!\n";
        #                                
        #                                system("rm $LastORFAln");
        #                        }
        #                }
        #        }else{
        #                print "The ORF $KnownORF is not present in $Qry\n";
        #                $NewReport -> [$j][$i+1] = "";
        #        }
        #        system("rm $ORFTemp");
        #}
        #open (FILE, ">$PresenceAbsence");
        #        for($j=0; $j<$LinesOnCoreGenome; $j++){
        #                for ($k=0; $k<$i+2; $k++){
        #                        if (defined ($NewReport -> [$j][$k])){
        #                                print FILE $NewReport -> [$j][$k], ",";      
        #                        }else{
        #                                $NewReport -> [$j][$k] = "";
        #                                print FILE $NewReport -> [$j][$k], 
        #                        }
        #                }
        #                print FILE "\n";
        #        }
        #close FILE;
#}

#@PresenceAbsenceArray = ReadFile($PresenceAbsence);
#($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, @PresenceAbsenceMatrix) = Matrix($PresenceAbsence);
#
#open (FILE,">$PresenceAbsence");
#        foreach $Line(@PresenceAbsenceArray){
#                @Fields = split(",",$Line);
#                chomp@Fields;
#                $Count = 0;
#                foreach $Element(@Fields){
#                        if ($Element ne ""){
#                                $Count++;
#                        }
#                }
#                if($Count == $ColumnsOnPresenceAbsence){
#                        print FILE "$Line\n";
#                }
#        }
#close FILE;
#
#system("mv $PresenceAbsence $CoreGenomeFile");
#($LinesOnCoreGenome, $ColumnsOnCoreGenome, @CoreGenomeMatrix) = Matrix($CoreGenomeFile);
#MakeDir($CoreSeqsPath);

for ($i=0; $i<$TotalQry; $i++){
       $Strain = $List[$i];
       $Db = $BlastPath ."/". $Strain;
       $OutCore = $CoreSeqsPath ."/". $Strain . "-CoreGenome" . $SeqExt;
       print "Building $Strain CoreGenome file...";
       for($j=1; $j<scalar@SharedORFs; $j++){ 
              $Gene = $Entry{$ORF}{$Strain};
              $cmd = `blastdbcmd -db $Db -dbtype nucl -entry "$Gene" -out $TempFile`;          
              $cmd = `cat $TempFile >> "$OutCore"`;   
              $cmd = `rm $TempFile`;
       }
       print "Done!\n";
}
exit;