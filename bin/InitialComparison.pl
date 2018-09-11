#!/usr/bin/perl -w
#################################################################################
#Scipt InitialComparison.pl                                                     #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          11 de abril de 2017                                             #
#################################################################################

use strict; 
use List::MoreUtils qw{any first_index};
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $List, $Trusted, $MolType, $eValue, $PIdent, $CPUs,
    $OutPath);

$Usage = "\tUsage: PreCoreGenome.pl <Project Name> <Genomes_List_File.ext> <TrustedORFeome.fasta> <Molecule type> <e-value> <Perc Ident> <CPUs> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName    = $ARGV[0];
$List           = $ARGV[1];
$Trusted        = $ARGV[2];
$MolType        = $ARGV[3];
$eValue         = $ARGV[4];
$PIdent         = $ARGV[5];
$CPUs           = $ARGV[6];
$OutPath        = $ARGV[7];

my ($Project, $ORFeomesPath, $ProjectGenomeList, $Sub, $Qry, $QryFile,
    $TrustedFile, $SubDb, $QryDb, $BlastReport, $cmd, $ORFsPath, $BlastPath,
    $TrustedORFeomePrefix, $SeqExt, $FastaAlnExt, $StoAlnExt, $HmmExt, $PanGenomeDb,
    $PresenceAbsence, $PanGenomeSeq, $Progress, $QryIDsFile, $TrustedIDsFile,
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
    $TrustedORFeomeDb, $QryExt, $ObservedPIdent, $Index, $LastCount, $Gene,
    $QryGenomeName, $QryGenomeRemainingORFs,$ORFSufix, $TestingORF, $ORFpath, $Hmm,
    $ORFTemp, $HmmSearch, $BestHit, $Entry, $NewCount, $PanGenomeDbDir,
    $PanGenomeHmmDb, $PurgeSearch, $ORFi
    );
my ($i, $j, $n, $k);
my (@List, @BlastReport, @ReportFields, @NonSharedQryIDs, @NonSharedTrustedIDs,
    @Duplicates, @QryIDs, @TrustedIDs, @SplitORFPath, @AnalizedORFs, @ID, @nHMMerReport,
    @BestHitArray, @Purge, @ORFs, @HmmScanHit);
my (%IDs);
my $OutReport = [ ];

$Project              = $OutPath ."/". $ProjectName;
$TrustedORFeomePrefix = Prefix($Trusted);

@List     = ReadFile($List);
$TotalQry = scalar@List;
$Qry      = $List[0];

$SeqExt = ".fasta";
$FastaAlnExt = ".aln.fasta";
$StoAlnExt = ".aln.sto";
$HmmExt = ".hmm";

if ($MolType eq "nucl"){
	$QryExt = ".ffn";
        $HmmSearch = "nhmmer";
}elsif($MolType eq "prot"){
	$QryExt = ".faa";
        $HmmSearch = "hmmsearch";
}

#Paths       
$ORFeomesPath           = $OutPath ."/". "ORFeomes";
$BlastPath              = $OutPath ."/". "Blast";
$ORFsPath               = $OutPath ."/". "ORFs";

#Inputs
$QryFile                = $ORFeomesPath ."/". $Qry . $QryExt;
$TrustedFile            = $ORFeomesPath ."/". $Trusted;
$QryDb                  = $BlastPath ."/". $Qry;
$PanGenomeDb            = $BlastPath ."/". "PanGenomeDb";
$TrustedORFeomeDb       = $BlastPath ."/". $TrustedORFeomePrefix . "Db";

#Outputs
$BlastReport            = $OutPath ."/". $ProjectName . "_UniqueBlastComparison.txt";
$PresenceAbsence        = $OutPath ."/". $ProjectName . "_Initial_Presence_Absence.csv";
$PanGenomeSeq           = $OutPath ."/". $ProjectName . "_PanGenome" . $SeqExt;
$Progress               = $OutPath ."/". $ProjectName . "_Progress.csv";
$Summary	        = $OutPath ."/". $ProjectName . "_Summary.txt";
$LogFile                = $OutPath ."/". $ProjectName . ".log";

MakeDir ($ORFsPath);
open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

@QryIDs = AnnotatedGenes($QryFile);
$TotalQryORFs = scalar@QryIDs;
@TrustedIDs = AnnotatedGenes($TrustedFile);
$TotalTrustedORFs = scalar@TrustedIDs;

$SharedQryORFsCounter = $TotalQryORFs;
$SharedTrustedORFsCounter = $TotalTrustedORFs;

print "\nLooking for the first shared ORFs between $TrustedORFeomePrefix and $Qry:\n";

if ($MolType eq "nucl"){
        $cmd = `blastn -query $QryFile -db $TrustedORFeomeDb -out $BlastReport -outfmt '6 qacc sacc length qlen slen qstart qend sstart send pident evalue bitscore' -evalue $eValue -max_hsps 1 -max_target_seqs 1 -qcov_hsp_perc 90 -perc_identity $PIdent -num_threads $CPUs`;
}elsif($MolType eq "prot"){
        $cmd = `blastp -query $QryFile -db $TrustedORFeomeDb -out $BlastReport -outfmt '6 qacc sacc length qlen slen qstart qend sstart send pident evalue bitscore' -evalue $eValue -max_hsps 1 -max_target_seqs 1 -qcov_hsp_perc 90 -num_threads $CPUs`;
}
	
@BlastReport = ReadFile($BlastReport);
$n = scalar@BlastReport;

$Count = 1;
$j = 0;
for ($i=0; $i<$n; $i++){
        $Counter = sprintf "%.4d", $Count;
        
        @ReportFields       = SplitTab($BlastReport[$i]);
	$QryId              = $ReportFields[0];
	$TrustedId          = $ReportFields[1];
        $ObservedPIdent     = $ReportFields[9];
        $QryOutFileName     = $Qry ."-". $QryId . $SeqExt;
        $TrustedOutFileName = $TrustedORFeomePrefix ."-". $TrustedId . $SeqExt;
        
        if ($PIdent < $ObservedPIdent){
        
                #GenesInBlastReport($QryIDsFile,$QryId); #Obtaining the names of those query genes that are shared with the trusted orfeome
                #GenesInBlastReport($TrustedIDsFile,$TrustedId); #Obtaining the names of those trusted genes that are shared with the query orfeome
                
                #@NonSharedQryIDs = DismissORFs($QryId, @NonSharedQryIDs); #Obtaining non shared query genes
                #@NonSharedTrustedIDs = DismissORFs($TrustedId, @NonSharedTrustedIDs); #Obtaining non shared trusted genes
                
                
                #Keep only non shared query genes
                if ( any { $_ eq $QryId} @QryIDs){  
                        $Index = first_index{$_ eq $QryId} @QryIDs;
                        splice@QryIDs,$Index,1;
                }
        
                #keep only non shared trusted genes
                if ( any { $_ eq $TrustedId} @TrustedIDs){  
                        $Index = first_index{$_ eq $TrustedId} @TrustedIDs;
                        splice@TrustedIDs,$Index,1;
                }
                
                $SharedORF = "ORF" ."_". $Counter;
                $SharedORFPath = $ORFsPath ."/". $SharedORF;
                
                $TrustedOut = $SharedORFPath ."/". $TrustedOutFileName;
                $QryOut = $SharedORFPath ."/". $QryOutFileName;
        
                $ToAlign = $SharedORFPath ."/". $SharedORF . $SeqExt;
                $fastaAln = $SharedORFPath ."/". "1-" . $SharedORF . $FastaAlnExt;
                $stoAln  = $SharedORFPath ."/". $SharedORF . $StoAlnExt;
                $OutHmm  = $SharedORFPath ."/". $SharedORF . $HmmExt;
                
                #$FindingTrusted = `find $ORFsPath -type f -name $TrustedOutFileName`;
                #if ($FindingTrusted ne ""){
                if (any {$_ eq $TrustedId} @AnalizedORFs){
                        
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
                        
                        push @AnalizedORFs, $TrustedId;
                
                        MakeDir($SharedORFPath);
                        Extract($TrustedORFeomePrefix,$TrustedORFeomeDb,$MolType,$TrustedId,$TrustedOut);
                        Extract($Qry,$QryDb,$MolType,$QryId,$QryOut);
                        Align($TrustedOut,$QryOut,$ToAlign,$fastaAln);
                        HMM($fastaAln,$MolType,$OutHmm,$CPUs);
                
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
}
$CoreGenomeSize = $j;

#Adding Non Shared ORFs and Pan-Genome building step
system("cp $TrustedFile $PanGenomeSeq");

print "\nProcessing Non Shared Genes\n"; #  <- Non shared genes from trusted orfeome
$TotalNonSharedORFs = scalar@TrustedIDs;
for($i=0; $i<$TotalNonSharedORFs; $i++){
        $NonSharedORFId = $TrustedIDs[$i];
        $NonSharedCounter = sprintf "%.4d",$CoreGenomeSize+$i+1;
        $NonSharedORF = "ORF" ."_". $NonSharedCounter;
        $NonSharedORFPath = $ORFsPath ."/". $NonSharedORF;
        $NonSharedORFSeq = $NonSharedORFPath ."/". $TrustedORFeomePrefix ."-". $NonSharedORFId . $SeqExt;
        $NonSharedORFHmm = $NonSharedORFPath ."/". $NonSharedORF . $HmmExt;
        $NonSharedORFAln = $NonSharedORFPath ."/". "1-" . $NonSharedORF . $FastaAlnExt;
        
        print "\nProcessing ORF $NonSharedCounter: \n";
        
        MakeDir($NonSharedORFPath);
        Extract($TrustedORFeomePrefix,$TrustedORFeomeDb,$MolType,$NonSharedORFId,$NonSharedORFSeq);
        $cmd = `cp $NonSharedORFSeq $NonSharedORFAln`;
        HMM($NonSharedORFAln,$MolType,$NonSharedORFHmm,$CPUs);
        
        $OutReport -> [$NonSharedCounter][0] = $NonSharedORF; 
        $OutReport -> [$NonSharedCounter][1] = $NonSharedORFId;
        $OutReport -> [$NonSharedCounter][2] = "";
        $Count++;
}

print "\nProcessing New ORFs\n"; # <- Non shared genes from query file
#Solving duplicates and fragmented genes
####################################
$TotalNewORFs = scalar@QryIDs;
@Purge = @QryIDs;
$PanGenomeDbDir = $OutPath ."/". "PanGenomeHmmDb";
MakeDir($PanGenomeDbDir);
for ($k=1; $k<$Count; $k++){
        $ORFSufix = sprintf "%.4d", $k;
        $TestingORF = "ORF" ."_". $ORFSufix;
        $ORFpath = $ORFsPath ."/". $TestingORF;
        $Hmm = $ORFpath ."/". $TestingORF . $HmmExt;
        $PanGenomeHmmDb = $PanGenomeDbDir ."/". "PanGenomeDb.hmm";
        system("cat $Hmm >> $PanGenomeHmmDb");
}
system ("hmmpress $PanGenomeHmmDb");
for($j=0;$j<$TotalNewORFs;$j++){
        $Gene = $QryIDs[$j];
        @ID = split("_",$Gene);
        $QryGenomeName = $ID[0];
        $QryGenomeRemainingORFs = $ORFeomesPath ."/". $Gene . "_Temp" . $SeqExt;
        $cmd = `blastdbcmd -db $QryDb -dbtype $MolType -entry "$Gene" > $QryGenomeRemainingORFs`;
        $PurgeSearch = $PanGenomeDbDir ."/". $Gene . "Purge.Temp";
        system ("hmmscan -E 1e-3 --cpu $CPUs --noali --tblout $PurgeSearch $PanGenomeHmmDb $QryGenomeRemainingORFs");
        @nHMMerReport = ReadFile($PurgeSearch);
        if(@nHMMerReport){
            if ( any { $_ eq $Gene} @Purge){
                $Index = first_index{$_ eq $Gene} @Purge;
                splice@Purge,$Index,1;
            }
        }
        system ("rm $PurgeSearch $QryGenomeRemainingORFs");
}
@QryIDs = @Purge;
####################################

$TotalNewORFs = scalar@QryIDs;
if ($TotalNewORFs > 0){
    for($i=0; $i<$TotalNewORFs; $i++){
            $NewORFId = $QryIDs[$i];
            $NewCounter = sprintf "%.4d",$NonSharedCounter+$i+1;
            $NewORF = "ORF" ."_". $NewCounter;
        	$NewORFPath = $ORFsPath ."/". $NewORF;
            $NewORFSeq = $NewORFPath ."/". $Qry ."-". $NewORFId . $SeqExt;
            $NewORFHmm = $NewORFPath ."/". $NewORF . $HmmExt;
            $NewORFAln = $NewORFPath ."/". "1-" . $NewORF . $FastaAlnExt;

            print "\nProcessing ORF $NewCounter: \n";

            MakeDir($NewORFPath);
            Extract($Qry,$QryDb,$MolType,$NewORFId,$NewORFSeq);
            $cmd = `cp $NewORFSeq $NewORFAln`;
            HMM($NewORFAln,$MolType,$NewORFHmm,$CPUs);
        
        	print "\tAdding ORF $NewCounter to PanGenome...";
            $cmd = `blastdbcmd -db $QryDb -dbtype $MolType -entry "$NewORFId" >> $PanGenomeSeq`;
        	print "Done!\n";
        
            $OutReport -> [$NewCounter][0] = $NewORF; 
            $OutReport -> [$NewCounter][1] = "";
            $OutReport -> [$NewCounter][2] = $NewORFId;
            $Count++;
    }
}

#Building a Presence/Absence genes file
open (FILE, ">>$PresenceAbsence");
for($i=0; $i<$NewCounter+1; $i++){
    for ($j=0; $j<3; $j++){
        print FILE $OutReport -> [$i][$j], ",";
    }
    print FILE "\n";
}
close FILE;

#Updating Trusted Orfeome Data Base
print "\nBuilding a Pan-Genome Data Base...";
$cmd = `makeblastdb -in $PanGenomeSeq -dbtype $MolType -parse_seqids -out $PanGenomeDb`;
print "Done!\n\n";

system ("rm -r $PanGenomeDbDir*");
MakeDir($PanGenomeDbDir);
for ($k=1; $k<$Count; $k++){
        $ORFSufix = sprintf "%.4d", $k;
        $TestingORF = "ORF" ."_". $ORFSufix;
        $ORFpath = $ORFsPath ."/". $TestingORF;
        $Hmm = $ORFpath ."/". $TestingORF . $HmmExt;
        $PanGenomeHmmDb = $PanGenomeDbDir ."/". "PanGenomeDb.hmm";
        system("cat $Hmm >> $PanGenomeHmmDb");
}

#Summary Report
$PanGenomeSize=`grep ">" $PanGenomeSeq | wc -l`;
chomp $PanGenomeSize;
open (FILE, ">$Summary");
        print FILE "Project Name: $ProjectName\n";
        print FILE "Included Genomes List File: $List\n";
        print FILE "Used CPU's: $CPUs\n";
        print FILE "e-value Threshold: $eValue\n";
        print FILE "Percentage of Identity: $PIdent\n";
        print FILE "Number of taxa: $TotalQry\n";        
        print FILE "Trusted ORFeome: $TotalTrustedORFs\n";
close FILE;

#Progress File
open (FILE, ">$Progress");
	print FILE "NumberOfNewStrains,AnalyzedStrain,PanGenome,CoreGenome,NewGenes\n";
	print FILE "0,$TrustedORFeomePrefix,$TotalTrustedORFs,$TotalTrustedORFs,0\n";
	print FILE "1,$Qry,$PanGenomeSize,$CoreGenomeSize,$TotalNewORFs\n";
close FILE;

exit;
