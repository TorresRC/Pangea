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
use lib '/Users/RC/lib';
use Routines;

my ($Usage, $ProjectName, $List, $TrustedORFeome, $CPUs);

$Usage = "\tUsage: PreCoreGenome.pl <Project Name> <Genomes_List_File.ext> <Path/To/TrustedORFeome.fasta> <CPUs>\n";
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
    $fastaAln, $SeqExt, $AlnExt, $HmmExt, $OutHmm, $PreCoreFile, $LogFile, $PanGenome,
    $DuplicateSubIDs, $QryIDsFile, $DuplicateBlastHits, $Duplicate, $stoExt, $stoAln,
    $PanGenomesDB, $TrustedORFeomePrefix, $QryIDs, $QryIdIndex, $NoTrustedGenes,
    $NoNewGenes, $c, $NewGene, $GeneId, $PanGenomeStats, $PanGenomeSize);
my (@List, @BlastReport, @ReportFields, @Duplicates, @QryIDs);
my (%IDs);
my $OutReport = [ ];

$MainPath               = "/Users/RC/CoreGenome/MacOS";
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
$BlastReport            = $Project ."/". $ProjectName ."-". $Qry ."_vs_". $TrustedORFeomePrefix ."_BlastReport.txt";
$PreCoreFile            = $Project ."/". $ProjectName . "_PreCoreGenome.csv";
$PanGenome              = $Project ."/". $ProjectName . "_PanGenome" . $SeqExt;
$PanGenomeStats         = $Project ."/". $ProjectName . "_PanGenomeStatistics.txt";
$QryIDsFile             = $Project ."/". $ProjectName . "_Shared" . $Qry . "GenesIDs.txt";
$DuplicateSubIDs        = $Project ."/". $ProjectName . "_ProbableDuplicateGenes.txt";
$DuplicateBlastHits     = $Project ."/". $ProjectName . "_DuplicateBlastHits.txt";
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
$n = scalar @BlastReport;
for ($i=0; $i<$n; $i++){
        $Counter = Counter($i);
        @ReportFields = SplitTab($BlastReport[$i]);
	$QryId = $ReportFields[0];
	$SubId = $ReportFields[1];
        
        open (FILE, ">>$QryIDsFile");
                print FILE "$QryId\n";
        close FILE;
        
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

	print "Analyzing ORF $Counter: \n";

	MakeDir($SharedORFPath); 

	print "\tExtracting ORF from $TrustedORFeomePrefix...";	
	$cmd = `blastdbcmd -db $PanGenomesDB -dbtype nucl -entry "$SubId" -out $SubOut`;
	print "Done!\n";

	print "\tExtracting ORF from $Qry...";
	$cmd = `blastdbcmd -db $QryDb -dbtype nucl -entry "$QryId" -out $QryOut`;
	print "Done!\n";

	print "\tAligning sequences...";
	$cmd = `cat $SubOut $QryOut > $ToAlign`;
	$cmd = `muscle -in $ToAlign -out $fastaAln -quiet`;
	print "Done!\n";

	print "\tBuilding a HMM...";
	$cmd = `hmmbuild --dna --cpu $CPUs $OutHmm $fastaAln`;
	print "Done!\n\n";
        
 	$OutReport -> [0][0] = "";		
        $OutReport -> [$i+1][0] = $SharedORF;
    
        $OutReport -> [0][1] = $TrustedORFeomePrefix;
        $OutReport -> [0][2] = $Qry;
        
        $OutReport -> [$i+1][1] = $SubId;
        $OutReport -> [$i+1][2] = $QryId;
}

#Duplicated Qry Shared Genes
system("uniq -d $QryIDsFile > $DuplicateSubIDs");

@Duplicates = ReadFile($DuplicateSubIDs);
chomp @Duplicates;
foreach $Duplicate(@Duplicates){
        system("grep -E $Duplicate $BlastReport >> $DuplicateBlastHits");
}

#Pan-Genome building step
print "\nBuilding Pan-Genome file\n";
system("cp $TrustedFile $PanGenome");
$NewCounter = $Counter+1;
$NoNewGenes = scalar@QryIDs;
for($c=0; $c<$NoNewGenes; $c++){
        $NewORFId = $QryIDs[$c];
        $NewORF = "ORF" ."_". $NewCounter+$a;
	$SharedORFPath = $ORFsPath ."/". $NewORF;
        MakeDir($SharedORFPath);  
        
        
        $cmd = `blastdbcmd -db $QryDb -dbtype nucl -entry "$NewORFId" >> $PanGenome`;
        Progress($NoNewGenes,$c);
}
$PanGenomeSize=`grep ">" $PanGenome | wc -l`;
chomp $PanGenomeSize;
open (FILE, ">$PanGenomeStats");
        print FILE "Old Pan-Genome: $NoTrustedGenes\nNew Genes: $NoNewGenes\nNew Pan-Genome: $PanGenomeSize";
close FILE;
print "Done!\n";

#PreCore-Genome building step
open (FILE, ">>$PreCoreFile");
for($a=0; $a<$n+1; $a++){
    for ($b=0; $b<3; $b++){
        print FILE $OutReport -> [$a][$b], ",";
    }
    print FILE "\n";
}
close FILE;

#Updating Trusted Orfeome Data Base
print "\nUpdating Trusted Orfeome Data Base with $NoNewGenes new genes...";
$cmd = `makeblastdb -in $PanGenome -dbtype nucl -parse_seqids -out $PanGenomesDB`;
print "Done!\n";

exit;