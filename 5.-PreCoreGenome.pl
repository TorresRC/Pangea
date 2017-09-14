#!/usr/bin/perl -w
#################################################################################
#   Programa Filter Features                                                    #
#   Nota: Se debe ajustar la ruta de lib y de la variable $PathSeq a las que    #
#   realmente se tengan donde se instalción el programa.                        #
#                                                                               #
# Programador:   Roberto C. Torres                                              #
# Fecha:  11 de abril de 2017                                                   #
#         6 de septiembre 2017 - Cambio de primer genoma por secuencia Trusted  #
#################################################################################

use strict; 
use List::MoreUtils qw{any};
use lib '/Users/rc/CoreGenome/src/lib';
use Routines;

my ($Usage, $ProjectName, $List, $TrustedORFeome, $eValue, $PIdent, $CPUs);

$Usage = "\tUsage: PreCoreGenome.pl <Project Name> <Genomes_List_File.ext> <TrustedORFeome.fasta> <e-value> <Perc Ident> <CPUs>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List = $ARGV[1];
$TrustedORFeome = $ARGV[2];
$eValue = $ARGV[3];
$PIdent = $ARGV[4];
$CPUs = $ARGV[5];

my ($MainPath, $Project, $ORFeomesPath, $ProjectGenomeList, $Sub, $Qry, $QryFile,
    $TrustedFile, $SubDb, $QryDb, $BlastReport, $cmd, $ORFsPath, $BlastPath, $Counter, $i,
    $n, $SubId, $QryId, $SharedORF, $SharedORFPath, $SubOut, $QryOut, $ToAlign, $ID,
    $fastaAln, $SeqExt, $AlnExt, $HmmExt, $OutHmm, $PanGenomeReport, $LogFile, $PanGenomeSeq,
    $DuplicatedTrustedIDs, $DuplicatedQryIDs, $QryIDsFile, $DuplicatedBlastHits, $Duplicate,
    $stoExt, $stoAln, $PanGenomesDB, $TrustedORFeomePrefix, $QryIDs, $QryIdIndex, $NoTrustedGenes,
    $NoNewGenes, $c, $NewGene, $GeneId, $Stats, $PanGenomeSize, $NewORFId, $NewCounter,
    $NewORF, $NewORFPath, $NewORFSeq, $NewORFAln, $NewORFHmm, $Summary,$CoreGenomeSize,
    $TrustedIDs, $TrustedIdIndex, $NoNonSharedORFs, $NonSharedORFId, $NonSharedCounter,
    $NonSharedORF, $NonSharedORFPath, $NonSharedORFSeq, $NonSharedORFHmm, $NonSharedORFAln, $x,
    $z, $d, $e, $TrustedIDsFile);
my (@List, @BlastReport, @ReportFields, @Duplicates, @QryIDs, @TrustedIDs);
my (%IDs);
my $OutReport = [ ];

$MainPath               = "/Users/rc/CoreGenome";
$Project                = $MainPath ."/". $ProjectName;
$ProjectGenomeList      = $Project ."/". $List;
@List = ReadFile($ProjectGenomeList);
$Qry  = $List[0];
$TrustedORFeomePrefix= Prefix($TrustedORFeome);

#File extentions 
$SeqExt = ".fasta";
$AlnExt = ".aln.fasta";
$stoExt = ".sto";
$HmmExt = ".hmm";

#Paths       
$ORFeomesPath           = $MainPath ."/". "ORFeomes";
$BlastPath              = $MainPath ."/". "Blast";
$ORFsPath               = $Project ."/". "ORFs";

#Inputs
$QryFile                = $ORFeomesPath ."/". $Qry . $SeqExt;
$TrustedFile            = $MainPath ."/". $TrustedORFeome;
$QryDb                  = $BlastPath ."/". $Qry;
$PanGenomesDB           = $BlastPath ."/". "PanGenomeDB";

#Output
$BlastReport            = $Project ."/". $ProjectName . "_UniqueBlastComparison.txt";
$PanGenomeReport        = $Project ."/". $ProjectName . "_Presence_Absence.csv";
$PanGenomeSeq           = $Project ."/". $ProjectName . "_PanGenome" . $SeqExt;
$Stats         		= $Project ."/". $ProjectName . "_Statistics.csv";
$QryIDsFile             = $Project ."/". $ProjectName . "_Shared_" . $Qry . "GenesIDs.txt";
$TrustedIDsFile         = $Project ."/". $ProjectName . "_Shared_" . $TrustedORFeomePrefix . "GenesIDs.txt";
$DuplicatedTrustedIDs   = $Project ."/". $ProjectName . "_BlastTrustedDuplicatedGenes.txt";
$DuplicatedQryIDs       = $Project ."/". $ProjectName . "_BlastQueryDuplicatedGenes.txt";
$DuplicatedBlastHits    = $Project ."/". $ProjectName . "_DuplicateBlastReport.txt";
$Summary		= $Project ."/". $ProjectName . "_Summary.txt";
$LogFile                = $Project ."/". $ProjectName . ".log";

#Numer of genes in the trusted file
$NoTrustedGenes = `grep ">" $TrustedFile | wc -l`;
chomp $NoTrustedGenes;

#Query gene IDs 
$QryIDs = `grep ">" $QryFile`;
$QryIDs =~ s/>//g;
@QryIDs = split('\n',$QryIDs);

#Trusted gene IDs
$TrustedIDs = `grep ">" $TrustedFile`;
$TrustedIDs =~ s/>//g;
@TrustedIDs = split('\n',$TrustedIDs);

MakeDir ($ORFsPath);
open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

print "\nLooking for the first shared ORFs between $TrustedORFeomePrefix and $Qry:\n";
$cmd = `blastn -query $QryFile -db $PanGenomesDB -out $BlastReport -outfmt '6 qacc sacc length qlen slen qstart qend sstart send pident evalue bitscore' -evalue $eValue -max_hsps 1 -max_target_seqs 1 -qcov_hsp_perc 90 -perc_identity $PIdent -num_threads $CPUs`;
	
@BlastReport = ReadFile($BlastReport);
$CoreGenomeSize = scalar@BlastReport;
$n = $CoreGenomeSize;
for ($i=0; $i<$n; $i++){
        $Counter = Counter($i);
        @ReportFields = SplitTab($BlastReport[$i]);
	$QryId = $ReportFields[0];
	$SubId = $ReportFields[1];
        
        
        #Obtaining Duplicate Blast Hits from both Query and Trusted sequences 
        open (FILE, ">>$QryIDsFile");
                print FILE "$QryId\n";
        close FILE;
        
        open (FILE, ">>$TrustedIDsFile");
                print FILE "$SubId\n";
        close FILE;
        
	#Keep only non query shared genes
        for($x=0;$x<scalar@QryIDs;$x++){
                $ID = $QryIDs[$x];
                $ID =~ s/\s//g;
                if($ID eq $QryId){
                        splice @QryIDs, $x, 1;
                }
        }

	#keep only non trusted shared genes
        for($z=0;$z<scalar@TrustedIDs; $z++){
                $ID = $TrustedIDs[$z];
                $ID =~ s/\s//g;
                if($ID eq $SubId){
                        splice @TrustedIDs, $z, 1;
                }
        }
        
	$SharedORF = "ORF" ."_". $Counter;
	$SharedORFPath = $ORFsPath ."/". $SharedORF;
	
	$SubOut = $SharedORFPath ."/". $TrustedORFeomePrefix ."-". $SubId . $SeqExt;
	$QryOut = $SharedORFPath ."/". $Qry ."-". $QryId . $SeqExt;

	$ToAlign = $SharedORFPath ."/". $SharedORF . $SeqExt;
	$fastaAln  = $SharedORFPath ."/". "1-" . $SharedORF . $AlnExt;
        $stoAln  = $SharedORFPath ."/". $SharedORF . $stoExt;
	$OutHmm  = $SharedORFPath ."/". $SharedORF . $HmmExt;

	print "\nAnalyzing ORF $Counter: \n";

	MakeDir($SharedORFPath);
        Extract($TrustedORFeomePrefix,$PanGenomesDB,$SubId,$SubOut);
        Extract($Qry,$QryDb,$QryId,$QryOut);
        Align($SubOut,$QryOut,$ToAlign,$fastaAln);
        HMM($CPUs,$OutHmm,$fastaAln);
        
 	$OutReport -> [0][0] = "";		
        $OutReport -> [$i+1][0] = $SharedORF;
    
        $OutReport -> [0][1] = $TrustedORFeomePrefix;
        $OutReport -> [0][2] = $Qry;
        
        $OutReport -> [$i+1][1] = $SubId;
        $OutReport -> [$i+1][2] = $QryId;
}

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
$NoNonSharedORFs = scalar@TrustedIDs;
for($a=0; $a<$NoNonSharedORFs; $a++){
        $NonSharedORFId = $TrustedIDs[$a];
        $NonSharedCounter = $Counter+$a+1;
        $NonSharedORF = "ORF" ."_". $NonSharedCounter;
        $NonSharedORFPath = $ORFsPath ."/". $NonSharedORF;
        $NonSharedORFSeq = $NonSharedORFPath ."/". $TrustedORFeomePrefix ."-". $NonSharedORFId . $SeqExt;
        $NonSharedORFHmm = $NonSharedORFPath ."/". $TrustedORFeomePrefix ."-". $NonSharedORFId . $HmmExt;
        $NonSharedORFAln = $NonSharedORFPath ."/". "1-" . $NonSharedORF . $AlnExt;
        
        print "\nProcessing ORF $NonSharedCounter: \n";
        
        MakeDir($NonSharedORFPath);
        Extract($TrustedORFeomePrefix,$PanGenomesDB,$NonSharedORFId,$NonSharedORFSeq);
        $cmd = `cp $NonSharedORFSeq $NonSharedORFAln`;
        HMM($CPUs,$NonSharedORFHmm,$NonSharedORFAln);
        
        $OutReport -> [$NonSharedCounter][0] = $NonSharedORF; 
        $OutReport -> [$NonSharedCounter][1] = $NonSharedORFId;
        $OutReport -> [$NonSharedCounter][2] = "";
}

print "\nProcessing New ORFs\n"; # <- Non shared genes from query file
$NoNewGenes = scalar@QryIDs;
for($b=0; $b<$NoNewGenes; $b++){
        $NewORFId = $QryIDs[$b];
        $NewCounter = $NonSharedCounter+$b+1;
        $NewORF = "ORF" ."_". $NewCounter;
	$NewORFPath = $ORFsPath ."/". $NewORF;
        $NewORFSeq = $NewORFPath ."/". $Qry ."-". $NewORFId . $SeqExt;
        $NewORFHmm = $NewORFPath ."/". $Qry ."-". $NewORFId . $HmmExt;
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

#Building a Presence amd Absence genes file
open (FILE, ">>$PanGenomeReport");
for($d=0; $d<$NewCounter+1; $d++){
    for ($e=0; $e<3; $e++){
        print FILE $OutReport -> [$d][$e], ",";
    }
    print FILE "\n";
}
close FILE;

#Updating Trusted Orfeome Data Base
print "\nUpdating Pan-Genome Data Base with $NoNewGenes new genes...";
$cmd = `makeblastdb -in $PanGenomeSeq -dbtype nucl -parse_seqids -out $PanGenomesDB`;
print "Done!\n\n";

#Summary Report
$PanGenomeSize=`grep ">" $PanGenomeSeq | wc -l`;
chomp $PanGenomeSize;
open (FILE, ">$Summary");
        print FILE "Trusted ORFeome: $NoTrustedGenes\nCore-Genome: $CoreGenomeSize\nPan-Genome: $PanGenomeSize\nNew Genes: $NoNewGenes";
close FILE;

#Statistics File
open (FILE, ">$Stats");
	print FILE "Number Of New Strains,Analyzed Strain,Core Genome,Pan Genome,New Genes\n";
	print FILE "0,$TrustedORFeomePrefix,$NoTrustedGenes,$NoTrustedGenes,0\n";
	print FILE "1,$Qry,$CoreGenomeSize,$PanGenomeSize,$NoNewGenes\n";
close FILE;

exit;

###############################################################################
sub Extract{
        my ($Qry, $DataBase,$Entry,$OutSeq, $null) = @_;
        print "\tExtracting ORF from $Qry...";	
        $cmd = `blastdbcmd -db $DataBase -dbtype nucl -entry "$Entry" -out $OutSeq`;
        print "Done!\n";
}

sub Align{
        my ($Seq1, $Seq2, $ToAlign, $AlnFile, $null) = @_;
	print "\tAligning sequences...";
	$cmd = `cat $Seq1 $Seq2 > $ToAlign`;
	$cmd = `muscle -in $ToAlign -out $AlnFile -quiet`;
	print "Done!\n";
}

sub HMM{
        my ($CPUs, $HmmFile, $AlnFile, $null) = @_;
	print "\tBuilding a HMM...";
	$cmd = `hmmbuild --dna --cpu $CPUs $HmmFile $AlnFile`;
	print "Done!\n";
}
