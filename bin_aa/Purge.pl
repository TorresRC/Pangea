#!/usr/bin/perl -w
#################################################################################
#Scipt ContigsLength.pl                                                         #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          11 de octubre de 2017                                           #
#################################################################################

use strict;
use Bio::SeqIO;
my %unique;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $List, $TrustedORFeome, $MainPath, $Sim);

$Usage = "\tUsage: FilterOrfeomes.pl <Project_Name> <List_File.ext> <Trusted_ORFeome.fasta> <Main_Path>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName    = $ARGV[0];
$List           = $ARGV[1];
$TrustedORFeome = $ARGV[2];
$MainPath       = $ARGV[3];
$Sim            = $ARGV[4];

my ($SeqExt, $Project, $MainList, $ORFeomesPath, $FilteredORFsPath, $TrustedFile,
    $FilteredTrustedORFeome, $LogFile, $Qry, $Unfiltered, $Filtered, $Redundant,
    $seqio, $outseq, $seqs, $id, $seq);
my ($n, $i);
my (@List);

$SeqExt                 = ".faa";
$Project                = $MainPath ."/". $ProjectName;
$MainList               = $Project ."/". $List;
$ORFeomesPath           = $Project ."/". "ORFeomes/Sorted";
$FilteredORFsPath       = $ORFeomesPath ."/". "Filtered";
$TrustedFile            = $ORFeomesPath ."/". $TrustedORFeome;
$FilteredTrustedORFeome = $FilteredORFsPath ."/". $TrustedORFeome;

$LogFile                = $Project ."/". $ProjectName . ".log";

#open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

MakeDir($FilteredORFsPath);


@List = ReadFile($MainList);
$n = scalar@List;

for ($i=0; $i<$n;$i++){
        $Qry = $List[$i] . $SeqExt;
        $Unfiltered = $ORFeomesPath ."/". $Qry;
        $Filtered = $ORFeomesPath ."/". $FilteredORFsPath ."/". $Qry;
        $Redundant = $ORFeomesPath ."/". $FilteredORFsPath ."/". $Qry . ".redundants";
        
        system("skipredundant -sprotein1 -mode 1 -threshold $Sim -gapopen 10.0 -gapextend 0.5 -outseq $Filtered -redundantoutseq $Redundant -auto");     
}


#$file   = $TrustedFile;
#$seqio  = Bio::SeqIO->new(-file => $file, -format => "fasta");
#$outseq = Bio::SeqIO->new(-file => ">$file.uniq", -format => "fasta");
#
#while($seqs = $seqio->next_seq) {
#  $id  = $seqs->display_id;
#  $seq = $seqs->seq;
#  unless(exists($unique{$seq})) {
#    $outseq->write_seq($seqs);
#    $unique{$seq} +=1;
#  }
#}

#for ($i=0; $i<$n; $i++){
#	$file   = $List[$i];
#        $seqio  = Bio::SeqIO->new(-file => $file, -format => "fasta");
#        $outseq = Bio::SeqIO->new(-file => ">$file.uniq", -format => "fasta");
#
#        while($seqs = $seqio->next_seq) {
#                $id  = $seqs->display_id;
#                $seq = $seqs->seq;
#                unless(exists($unique{$seq})) {
#                        $outseq->write_seq($seqs);
#                        $unique{$seq} +=1;
#                }
#        }
#}
