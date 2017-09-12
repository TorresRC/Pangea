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
use lib '/home/rtorres/CoreGenome/src/lib';
use Routines;

my ($Usage, $ProjectName, $List, $TrustedORFeome, $CPUs);

$Usage = "\tUsage: PreCoreGenome.pl <Project Name> <Genomes_List_File.ext> <TrustedORFeome.fasta> <CPUs>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List = $ARGV[1];
$TrustedORFeome = $ARGV[2];
$CPUs = $ARGV[3];

my ($MainPath, $Project, $ORFeomesPath, $ProjectGenomeList, $Sub, $Qry, $QryFile,
    $TrustedFile, $SubDb, $QryDb, $BlastReport, $cmd, $ORFsPath, $BlastPath, $Counter, $i,
    $n, $SubId, $QryId, $SharedORF, $SharedORFPath, $SubOut, $QryOut, $ToAlign, $ID,
    $fastaAln, $SeqExt, $AlnExt, $HmmExt, $OutHmm, $PanGenomeReport, $LogFile, $PanGenomeSeq,
    $DuplicateSubIDs, $QryIDsFile, $DuplicateBlastHits, $Duplicate, $stoExt, $stoAln,
    $PanGenomesDB, $TrustedORFeomePrefix, $QryIDs, $QryIdIndex, $NoTrustedGenes,
    $NoNewGenes, $c, $NewGene, $GeneId, $Stats, $PanGenomeSize, $NewORFId, $NewCounter,
    $NewORF, $NewORFPath, $NewORFSeq, $NewORFAln, $NewORFHmm, $Summary,$CoreGenomeSize);
my (@List, @BlastReport, @ReportFields, @Duplicates, @QryIDs);
my (%IDs);
my $OutReport = [ ];

$MainPath               = "/home/rtorres/CoreGenome";
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
#$BlastReport           = $Project ."/". $ProjectName ."-". $Qry ."_vs_". $TrustedORFeomePrefix ."_BlastReport.txt";
$BlastReport            = $Project ."/". $ProjectName . "_UniqueBlastComparison.txt";
$PanGenomeReport        = $Project ."/". $ProjectName . "_PanGenome.csv";
$PanGenomeSeq           = $Project ."/". $ProjectName . "_PanGenome" . $SeqExt;
$Stats         		= $Project ."/". $ProjectName . "_Statistics.csv";
$QryIDsFile             = $Project ."/". $ProjectName . "_Shared_" . $Qry . "GenesIDs.txt";
$DuplicateSubIDs        = $Project ."/". $ProjectName . "_ProbableDuplicateGenes.txt";
$DuplicateBlastHits     = $Project ."/". $ProjectName . "_DuplicateBlastHits.txt";
$Summary		= $Project ."/". $ProjectName . "_Summary.txt";
$LogFile                = $Project ."/". $ProjectName . ".log";

#Numer of genes in the trusted file
$NoTrustedGenes = `grep ">" $TrustedFile | wc -l`;
chomp $NoTrustedGenes;

#Complete query genes ID 
$QryIDs = `grep ">" $QryFile`;
$QryIDs =~ s/>//g;
@QryIDs = split('\n',$QryIDs);

MakeDir ($ORFsPath);
open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

print "\nLooking for the first shared ORFs between $TrustedORFeomePrefix and $Qry:\n";
$cmd = `blastn -query $QryFile -db $PanGenomesDB -out $BlastReport -outfmt '6 qacc sacc qlen slen length qstart qend sstart send pident evalue bitscore' -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 -qcov_hsp_perc 90 -perc_identity 80 -num_threads $CPUs`;
	
@BlastReport = ReadFile($BlastReport);
$CoreGenomeSize = scalar@BlastReport;
$n = $CoreGenomeSize;
for ($i=0; $i<$n; $i++){
        $Counter = Counter($i);
        @ReportFields = SplitTab($BlastReport[$i]);
	$QryId = $ReportFields[0];
	$SubId = $ReportFields[1];
        
        open (FILE, ">>$QryIDsFile");
                print FILE "$QryId\n";
        close FILE;
        
	#Removing Shared Genes from Query
        $QryIdIndex = grep {$QryIDs[$_] eq $QryId} 0..$#QryIDs;
        if (any {$_ eq $QryId} @QryIDs){
                splice @QryIDs, $QryIdIndex-1, 1;
        }

	#Removing Shared Genes from Trusted ORFeome
	$QryIdIndex = grep {$QryIDs[$_] eq $QryId} 0..$#QryIDs;
        if (any {$_ eq $QryId} @QryIDs){
                splice @QryIDs, $QryIdIndex-1, 1;
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

#Duplicated Qry Genes
system("uniq -d $QryIDsFile > $DuplicateSubIDs");

@Duplicates = ReadFile($DuplicateSubIDs);
chomp @Duplicates;
foreach $Duplicate(@Duplicates){
        system("grep -E $Duplicate $BlastReport >> $DuplicateBlastHits");
}

#Adding Non Shared ORFs and Pan-Genome building step
print "\nAnalyzing New ORFs\n";
system("cp $TrustedFile $PanGenomeSeq");
$NoNewGenes = scalar@QryIDs;
for($a=0; $a<$NoNewGenes; $a++){
        $NewORFId = $QryIDs[$a];
        $NewCounter = $Counter+$a+1;
        $NewORF = "ORF" ."_". $NewCounter;
	$NewORFPath = $ORFsPath ."/". $NewORF;
        $NewORFSeq = $NewORFPath ."/". $Qry ."-". $NewORFId . $SeqExt;
        $NewORFHmm = $NewORFPath ."/". $Qry ."-". $NewORFId . $HmmExt;
        $NewORFAln = $NewORFPath ."/". "1-" . $NewORF . $AlnExt;

        print "\nAnalyzing ORF $NewCounter: \n";

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

#Building a Pan-Genome Report
open (FILE, ">>$PanGenomeReport");
for($b=0; $b<$NewCounter+1; $b++){
    for ($c=0; $c<3; $c++){
        print FILE $OutReport -> [$b][$c], ",";
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

exit;
