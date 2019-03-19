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

my ($Usage, $ProjectName, $Columns, $CoreGenome, $ORFsPath, $OutPath);

$Usage = "\tUsage: ConsensusPanGenome.pl <Project Name> <Columns> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$Columns     = $ARGV[1];
$CoreGenome  = $ARGV[2];
$ORFsPath    = $ARGV[3];
$OutPath     = $ARGV[4];

my($Line, $ORF, $ORFAln, $Name, $Seq, $Key, $AlnExt, $ORFPath, $Genome,
   $AlignedCoresPath, $SeqExt, $AlignedCore, $CoreFile, $Header);
my($i, $j, $LinesOnCoreGenome, $ColumnsOnCoreGenome);
my(@CoreGenome, @CoreGenomeFields, @CoreGenomeArray, @File, @IndexedName, @Seq);
my(%Seq);

$AlignedCoresPath = $OutPath ."/". "CoreSequences";
$SeqExt = ".fasta";
$AlnExt = ".aln" . $SeqExt;

#print "Loading the Core-Genome table:\n";
($LinesOnCoreGenome, $ColumnsOnCoreGenome, @CoreGenomeArray) = Matrix($CoreGenome);

print "Buiding aligned Core-Genomes fasta sequences:\n";
for ($i=1; $i<$ColumnsOnCoreGenome; $i++){
	$Genome = $CoreGenomeArray[0]->[$i];
	$AlignedCore = $AlignedCoresPath ."/". $Genome . "-CoreGenome". $AlnExt;
	open (FILE, ">>$AlignedCore");
		print FILE ">$Genome\n";
	close FILE;
}
for ($i=1; $i<$LinesOnCoreGenome; $i++){
        $ORF = $CoreGenomeArray[$i]->[0];
        $ORFPath = $ORFsPath ."/". $ORF;
        $ORFAln = `find $ORFPath -type f -name \*$AlnExt`;
        chomp $ORFAln;
        %Seq = ();
        
        open (FILE, $ORFAln) or die ("Could not open $ORFAln file. On $0 line 52.");
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
                        #print FILE "X\n";
			print FILE "\n";
                close FILE;
        }
        Progress($LinesOnCoreGenome, $i);
}

#print "Refining Core-Genome files:\n";
#for ($i=2; $i<$ColumnsOnCoreGenome; $i++){
#        $Header = $CoreGenomeArray[0]->[$i];
#        $CoreFile = $AlignedCoresPath ."/". $Header . "-CoreGenome". $AlnExt;
#        @Seq = ReadFile($CoreFile);
#        open (FILE, ">$CoreFile");
#                print FILE ">$Header\n";
#                foreach $Line (@Seq){
#                        print FILE "$Line\n";
#                }
#        close FILE;
#        Progress($ColumnsOnCoreGenome, $i);
#}
exit;
