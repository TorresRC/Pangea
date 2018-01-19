#!/usr/bin/perl -w
#################################################################################
#Scipt ContigsLength.pl                                                         #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          11 de abril de 2017                                             #
#################################################################################

use strict; 
use List::MoreUtils qw{any};
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $List, $TrustedORFeome, $eValue, $PIdent, $CPUs, $MainPath);

$Usage = "\tUsage: PreCoreGenome.pl <Project Name> <Genomes_List_File.ext> <TrustedORFeome.fasta> <e-value> <Perc Ident> <CPUs> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName    = $ARGV[0];
$List           = $ARGV[1];
$TrustedORFeome = $ARGV[2];
$eValue         = $ARGV[3];
$PIdent         = $ARGV[4];
$CPUs           = $ARGV[5];
$MainPath       = $ARGV[6];

my ($Project, $ORFeomesPath, $ProjectGenomeList, $Sub, $Qry, $QryFile,
    $TrustedFile, $SubDb, $QryDb, $BlastReport, $cmd, $ORFsPath, $BlastPath,
    $TrustedORFeomePrefix, $SeqExt, $AlnExt, $stoExt, $HmmExt, $PanGenomeDb,
    $PresenceAbsence, $PanGenomeSeq, $Stats, $QryIDsFile, $TrustedIDsFile,
    $DuplicatedTrustedIDs, $DuplicatedQryIDs, $DuplicatedBlastHits, $Summary,
    $LogFile, $TotalQryORFs, $TotalTrustedORFs, $CoreGenomeSize, $Count,
    $Counter, $QryId, $TrustedId, $ID, $SharedORF, $SharedORFPath,
    $TrustedOutFileName, $TrustedOut, $QryOutFileName, $QryOut, $ToAlign,
    $fastaAln, $stoAln, $OutHmm, $FindingTrusted, $KnownORFIndex, $KnownORFName,
    $KnownORFPath, $FileToKnownORF, $KnownFastaFile, $TempKnownFastaFile,
    $KnownAln, $KnownHmm, $Duplicate, $TotalNonSharedORFs, $NonSharedORFId,
    $NonSharedCounter, $NonSharedORF, $NonSharedORFPath, $NonSharedORFHmm,
    $NonSharedORFSeq, $NonSharedORFAln, $TotalNewORFs, $NewORFId, $NewCounter,
    $NewORF, $NewORFPath, $NewORFSeq, $NewORFHmm, $NewORFAln, $PanGenomeSize,
    $SharedQryORFsCounter, $SharedTrustedORFsCounter, $TotalQry,
    $TrustedORFeomeDb);
my ($a, $d, $e, $i, $j, $n, $x, $z);
my (@List, @BlastReport, @ReportFields, @NonSharedQryIDs, @NonSharedTrustedIDs,
    @Duplicates, @QryIDs, @TrustedIDs, @SplitORFPath);
my (%IDs);
my $OutReport = [ ];

$Project                = $MainPath ."/". $ProjectName;
$ProjectGenomeList      = $Project ."/". $List;
@List = ReadFile($ProjectGenomeList);
$TotalQry = scalar@List;
$Qry  = $List[0];
$TrustedORFeomePrefix= Prefix($TrustedORFeome);

#File extentions 
$SeqExt = ".fasta";
$AlnExt = ".aln.fasta";
$stoExt = ".sto";
$HmmExt = ".hmm";

#Paths       
$ORFeomesPath           = $Project ."/". "ORFeomes" ."/". "Sorted" ."/". "Filtered";
$BlastPath              = $Project ."/". "Blast";
$ORFsPath               = $Project ."/". "ORFs";

#Inputs
$QryFile                = $ORFeomesPath ."/". $Qry . ".ffn";
$TrustedFile            = $ORFeomesPath ."/". $TrustedORFeome;
$QryDb                  = $BlastPath ."/". $Qry;
$PanGenomeDb            = $BlastPath ."/". "PanGenomeDb";
$TrustedORFeomeDb       = $BlastPath ."/". $TrustedORFeomePrefix . "Db";

#Output
$BlastReport            = $Project ."/". $ProjectName . "_UniqueBlastComparison.txt";
$PresenceAbsence        = $Project ."/". $ProjectName . "_Initial_Presence_Absence.csv";
$PanGenomeSeq           = $Project ."/". $ProjectName . "_PanGenome" . $SeqExt;
$Stats         		    = $Project ."/". $ProjectName . "_Statistics.csv";
$QryIDsFile             = $Project ."/". $ProjectName . "_Shared_" . $Qry . "GenesIDs.txt";
$TrustedIDsFile         = $Project ."/". $ProjectName . "_Shared_" . $TrustedORFeomePrefix . "GenesIDs.txt";
$DuplicatedTrustedIDs   = $Project ."/". $ProjectName . "_BlastTrustedDuplicatedGenes.txt";
$DuplicatedQryIDs       = $Project ."/". $ProjectName . "_BlastQueryDuplicatedGenes.txt";
$DuplicatedBlastHits    = $Project ."/". $ProjectName . "_DuplicateBlastReport.txt";
$Summary	        	= $Project ."/". $ProjectName . "_Summary.txt";
$LogFile                = $Project ."/". $ProjectName . ".log";

MakeDir ($ORFsPath);
open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

@QryIDs = AnnotatedGenes($QryFile); #Get all names of annotated genes from query  
@TrustedIDs = AnnotatedGenes($TrustedFile); #Get all names of annotated genes from the trusted orfeome
$TotalQryORFs = scalar@QryIDs;
$TotalTrustedORFs = scalar@TrustedIDs;
#@NonSharedQryIDs = @QryIDs;
#@NonSharedTrustedIDs = @TrustedIDs;
$SharedQryORFsCounter = $TotalQryORFs;
$SharedTrustedORFsCounter = $TotalTrustedORFs;

print "\nLooking for the first shared ORFs between $TrustedORFeomePrefix and $Qry:\n";
$cmd = `blastn -query $QryFile -db $TrustedORFeomeDb -out $BlastReport -outfmt '6 qacc sacc length qlen slen qstart qend sstart send pident evalue bitscore' -evalue $eValue -max_hsps 1 -max_target_seqs 1 -qcov_hsp_perc 90 -perc_identity $PIdent -num_threads $CPUs`;
	
@BlastReport = ReadFile($BlastReport);
$n = scalar@BlastReport;
$Count = 1;
$j = 0;
for ($i=0; $i<$n; $i++){
        $Counter = sprintf "%.4d", $Count;
        
        @ReportFields = SplitTab($BlastReport[$i]);
	$QryId = $ReportFields[0]; #Id of the query gene in the blast hit
	$TrustedId = $ReportFields[1]; #Id of the trusted gene in the blast hit
        $QryOutFileName = $Qry ."-". $QryId . $SeqExt;
        $TrustedOutFileName = $TrustedORFeomePrefix ."-". $TrustedId . $SeqExt;
        
        GenesInBlastReport($QryIDsFile,$QryId); #Obtaining the names of those query genes that are shared with the trusted orfeome
        GenesInBlastReport($TrustedIDsFile,$TrustedId); #Obtaining the names of those trusted genes that are shared with the query orfeome
        
        #@NonSharedQryIDs = DismissORFs($QryId, @NonSharedQryIDs); #Obtaining non shared query genes
        #@NonSharedTrustedIDs = DismissORFs($TrustedId, @NonSharedTrustedIDs); #Obtaining non shared trusted genes
        
        #Keep only non shared query genes
        for($x=0;$x<$SharedQryORFsCounter;$x++){
                $ID = $QryIDs[$x];
                if($ID eq $QryId){
                        splice @QryIDs, $x, 1;
                        $SharedQryORFsCounter--;
                }
        }

	#keep only non shared trusted genes
        for($z=0;$z<$SharedTrustedORFsCounter; $z++){
                $ID = $TrustedIDs[$z];
                if($ID eq $TrustedId){
                        splice @TrustedIDs, $z, 1;
                        $SharedTrustedORFsCounter--;
                }
        }
        
	$SharedORF = "ORF" ."_". $Counter;
	$SharedORFPath = $ORFsPath ."/". $SharedORF;
	
        $TrustedOut = $SharedORFPath ."/". $TrustedOutFileName;
        $QryOut = $SharedORFPath ."/". $QryOutFileName;

	$ToAlign = $SharedORFPath ."/". $SharedORF . $SeqExt;
	$fastaAln  = $SharedORFPath ."/". "1-" . $SharedORF . $AlnExt;
        $stoAln  = $SharedORFPath ."/". $SharedORF . $stoExt;
	$OutHmm  = $SharedORFPath ."/". $SharedORF . $HmmExt;
        
        $FindingTrusted = `find $ORFsPath -type f -name $TrustedOutFileName`;
        if ($FindingTrusted ne ""){
                
                #@SplitORFPath = split("/",$FindingTrusted);
                #$KnownORFIndex = $#SplitORFPath-1;
                #$KnownORFName = $SplitORFPath[$KnownORFIndex];
                #$KnownORFPath = $ORFsPath ."/". $KnownORFName;
                #$FileToKnownORF = $KnownORFPath ."/". $QryOutFileName;
                #$KnownFastaFile = $KnownORFPath ."/". $KnownORFName . $SeqExt;
                #$TempKnownFastaFile = $KnownFastaFile . ".temp";
                #$cmd = `cp $KnownFastaFile $TempKnownFastaFile`;
                #$KnownAln  = $KnownORFPath ."/". "1-" . $KnownORFName . $AlnExt;
                #$KnownHmm = $KnownORFPath ."/". $KnownORFName . $HmmExt;
                
                #Extract($Qry,$QryDb,$QryId,$FileToKnownORF);
                #Align($TempKnownFastaFile,$FileToKnownORF,$KnownFastaFile,$KnownAln);
                #HMM($CPUs,$KnownHmm,$KnownAln);

                #$cmd = `rm $TempKnownFastaFile`;
        }else{
        	print "\nAnalyzing ORF $Counter: \n";
        
        	MakeDir($SharedORFPath);
                Extract($TrustedORFeomePrefix,$TrustedORFeomeDb,$TrustedId,$TrustedOut);
                Extract($Qry,$QryDb,$QryId,$QryOut);
                Align($TrustedOut,$QryOut,$ToAlign,$fastaAln);
                HMM($CPUs,$OutHmm,$fastaAln);
        
                $OutReport -> [0][0] = "ORF";		
                $OutReport -> [$j+1][0] = $SharedORF;
                
                $OutReport -> [0][1] = $TrustedORFeomePrefix;
                $OutReport -> [0][2] = $Qry;
                
                $OutReport -> [$j+1][1] = $TrustedId;
                $OutReport -> [$j+1][2] = $QryId;
                
                $Count++;
                $j++;
        }
}
$CoreGenomeSize = $j;

#Obtaining Duplicated Genes Files 
system("uniq -d $TrustedIDsFile > $DuplicatedTrustedIDs");
@Duplicates = ReadFile($DuplicatedTrustedIDs);
chomp @Duplicates;
foreach $Duplicate(@Duplicates){
        system("grep -E $Duplicate $BlastReport >> $DuplicatedBlastHits");
}

system("uniq -d $QryIDsFile > $DuplicatedQryIDs");
@Duplicates = ReadFile($DuplicatedQryIDs);
chomp @Duplicates;
foreach $Duplicate(@Duplicates){
        system("grep -E $Duplicate $BlastReport >> $DuplicatedBlastHits");
}

#Adding Non Shared ORFs and Pan-Genome building step
system("cp $TrustedFile $PanGenomeSeq");

print "\nProcessing Non Shared Genes\n"; #  <- Non shared genes from trusted orfeome
#$TotalNonSharedORFs = scalar@NonSharedTrustedIDs;
$TotalNonSharedORFs = scalar@TrustedIDs;
for($a=0; $a<$TotalNonSharedORFs; $a++){
        #$NonSharedORFId = $NonSharedTrustedIDs[$a];
        $NonSharedORFId = $TrustedIDs[$a];
        $NonSharedCounter = sprintf "%.4d",$CoreGenomeSize+$a+1;
        $NonSharedORF = "ORF" ."_". $NonSharedCounter;
        $NonSharedORFPath = $ORFsPath ."/". $NonSharedORF;
        $NonSharedORFSeq = $NonSharedORFPath ."/". $TrustedORFeomePrefix ."-". $NonSharedORFId . $SeqExt;
        $NonSharedORFHmm = $NonSharedORFPath ."/". $NonSharedORF . $HmmExt;
        $NonSharedORFAln = $NonSharedORFPath ."/". "1-" . $NonSharedORF . $AlnExt;
        
        print "\nProcessing ORF $NonSharedCounter: \n";
        
        MakeDir($NonSharedORFPath);
        Extract($TrustedORFeomePrefix,$TrustedORFeomeDb,$NonSharedORFId,$NonSharedORFSeq);
        $cmd = `cp $NonSharedORFSeq $NonSharedORFAln`;
        HMM($CPUs,$NonSharedORFHmm,$NonSharedORFAln);
        
        $OutReport -> [$NonSharedCounter][0] = $NonSharedORF; 
        $OutReport -> [$NonSharedCounter][1] = $NonSharedORFId;
        $OutReport -> [$NonSharedCounter][2] = "";
}

print "\nProcessing New ORFs\n"; # <- Non shared genes from query file
#$TotalNewORFs = scalar@NonSharedQryIDs;
$TotalNewORFs = scalar@QryIDs;
for($b=0; $b<$TotalNewORFs; $b++){
        #$NewORFId = $NonSharedQryIDs[$b];
        $NewORFId = $QryIDs[$b];
        $NewCounter = sprintf "%.4d",$NonSharedCounter+$b+1;
        $NewORF = "ORF" ."_". $NewCounter;
	$NewORFPath = $ORFsPath ."/". $NewORF;
        $NewORFSeq = $NewORFPath ."/". $Qry ."-". $NewORFId . $SeqExt;
        $NewORFHmm = $NewORFPath ."/". $NewORF . $HmmExt;
        $NewORFAln = $NewORFPath ."/". "1-" . $NewORF . $AlnExt;

        print "\nProcessing ORF $NewCounter: \n";

        MakeDir($NewORFPath);
        Extract($Qry,$QryDb,$NewORFId,$NewORFSeq);
        $cmd = `cp $NewORFSeq $NewORFAln`;
        HMM($CPUs,$NewORFHmm,$NewORFAln);
        
	print "\tAdding ORF $NewCounter to PanGenome...";
        $cmd = `blastdbcmd -db $QryDb -dbtype nucl -entry "$NewORFId" >> $PanGenomeSeq`;
	print "Done!\n";
        
        $OutReport -> [$NewCounter][0] = $NewORF; 
        $OutReport -> [$NewCounter][1] = "";
        $OutReport -> [$NewCounter][2] = $NewORFId;
}

#Building a Presence/Absence genes file
open (FILE, ">>$PresenceAbsence");
for($d=0; $d<$NewCounter+1; $d++){
    for ($e=0; $e<3; $e++){
        print FILE $OutReport -> [$d][$e], ",";
    }
    print FILE "\n";
}
close FILE;

#Updating Trusted Orfeome Data Base
print "\nBuilding a Pan-Genome Data Base...";
$cmd = `makeblastdb -in $PanGenomeSeq -dbtype nucl -parse_seqids -out $PanGenomeDb`;
print "Done!\n\n";

#Summary Report
$PanGenomeSize=`grep ">" $PanGenomeSeq | wc -l`;
chomp $PanGenomeSize;
open (FILE, ">$Summary");
        print FILE "Project Name: $ProjectName\n";
        print FILE "Included Genomes List File: $ProjectGenomeList\n";
        print FILE "Used CPU's: $CPUs\n";
        print FILE "e-value Threshold: $eValue\n";
        print FILE "Percentage of Identity: $PIdent\n";
        print FILE "Number of taxa: $TotalQry\n";        
        print FILE "Trusted ORFeome: $TotalTrustedORFs\n";
close FILE;

#Statistics File
open (FILE, ">$Stats");
	print FILE "NumberOfNewStrains,AnalyzedStrain,CoreGenome,PanGenome,NewGenes\n";
	print FILE "0,$TrustedORFeomePrefix,$TotalTrustedORFs,$TotalTrustedORFs,0\n";
	print FILE "1,$Qry,$CoreGenomeSize,$PanGenomeSize,$TotalNewORFs\n";
close FILE;

exit;
