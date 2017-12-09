#!/usr/bin/perl -w
#################################################################################
#Scipt CoreAlign.pl                                                             #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          19 de octubre de 2017                                           #
#################################################################################
use strict; 
use List::MoreUtils qw{any};
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $Columns, $MainPath);

$Usage = "\tUsage: ConsensusPanGenome.pl <Project Name> <Columns> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$Columns     = $ARGV[1];
$MainPath    = $ARGV[2];

my($Project, $ORFsPath, $CoreGenome, $Line, $ORF, $ORFAln, $Name, $Seq, $Key,
   $AlnExt, $ORFPath, $Genome, $AlignedCoresPath, $SeqExt, $AlignedCore);
my($i, $j, $LinesOnCoreGenome, $ColumnsOnCoreGenome);
my(@CoreGenome, @CoreGenomeFields, @CoreGenomeArray, @File, @IndexedName);
my(%Seq);

$Project = $MainPath ."/". $ProjectName;
$ORFsPath = $Project ."/". "ORFs";
$CoreGenome = $Project ."/". $ProjectName . "_CoreGenome.csv";
$AlignedCoresPath = $Project ."/". "CoreSequences";
$SeqExt = ".fasta";
$AlnExt = ".aln" . $SeqExt;

@CoreGenome = ReadFile($CoreGenome);
$LinesOnCoreGenome = scalar@CoreGenome;
 
print "Loading the Core-Genome table:\n";
for ($i=0; $i<$LinesOnCoreGenome; $i++){
     $Line = $CoreGenome[$i];
     @CoreGenomeFields = split(",",$Line);
     $ColumnsOnCoreGenome = scalar@CoreGenomeFields;
     push (@CoreGenomeArray, [@CoreGenomeFields]);
     Progress($LinesOnCoreGenome, $i);
}

print "Buiding aligned Core-Genomes fasta sequences:\n";
for ($i=1; $i<$LinesOnCoreGenome; $i++){
        $ORF = $CoreGenomeArray[$i]->[0];
        $ORFPath = $ORFsPath ."/". $ORF;
        $ORFAln = `find $ORFPath -type f -name \*$AlnExt`;
        chomp $ORFAln;
        %Seq = ();
        
        #print "\nAligning $ORF:\n";
        
        open (FILE, $ORFAln) or die ("Could not open $ORFAln file. On $0 line 63.");
        while (<FILE>){
                next unless /\S/;
                next if /^\s*\#/;
                if (/^\s*\/\//) {
                }else{
                        chomp;
                        ($Name, $Seq) = split;
                        $Seq =~ s/\./-/g;
                        $Seq{$Name} .= $Seq;
                }
        }
        close FILE;
        
        foreach $Key (sort keys %Seq){
                @IndexedName = split("_",$Key);
                $Genome = $IndexedName[0];
                $AlignedCore = $AlignedCoresPath ."/". $Genome . "-CoreGenome". $AlnExt;
                open (FILE, ">>$AlignedCore");
                        #print FILE ">$ORF~$Key\n";
                        for ($j=0; $j<length$Seq{$Key}; $j+=$Columns){
                                print FILE substr($Seq{$Key}, $j, $Columns), "\n";
                        }
                        print FILE "X\n";
                close FILE;
        }
        Progress($LinesOnCoreGenome, $i);
}
exit;
