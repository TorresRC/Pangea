#!/usr/bin/perl -w
#################################################################################
#                                                                               #
# Author: Roberto C. Torres, PhD (torres.roberto.c@gmail.com  )                 #
# Date: September 17th, 2019                                                    #
#                                                                               #
#################################################################################
use strict;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my $Src      = "$FindBin::Bin";

my $Usage    = "USAGE:\n  $0 [--help] [--vcf file] [--snpaln file] [--out path]
\n  Use \'--help\' to print detailed descriptions of options.\n\n";

my ($Help, $PresenceAbsence, $GenomesPath, $ORFsPath, $eValue, $OutPath);

GetOptions(
        'help'         => \$Help,
        'presabs|m=s'  => \$PresenceAbsence,
        'genomes|g=s'  => \$GenomesPath,
        'orfs|r=s'     => \$ORFsPath,
        'evalue|e=s'   => \$eValue,
        'out|o=s'      => \$OutPath
        ) or die $Usage;

my ($ORF, $Sample, $PanOrfHmm, $PanOrfAln, $SampleGenome, $OrfOutDir,
    $TempOutDir, $HmmHits, $SampleGeneFile, $CoreOrfAln, $CoreOrfHmm, $OutDir,
    $SingleSeq, $Header, $Seq, $Hmm, $Aln, $BestHit, $ContigName, $Coord1,
    $Coord2, $GeneStart, $GeneEnd, $GeneLen, $GeneSeq, $OutSeq, $SampleOrfAln,
    $TempAln);
my ($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence,$i, $j);
my (@PresenceAbsenceMatrix, @HmmHits, @BestHitData, @MultifastaSampleGenome);
my (%Seq, %AlnSeq, %seq, %CatSeq);

my $CoordsReport = [ ];
$CoordsReport -> [0][0] = "";

my $Sto2FastaScript = $Src . "/stockholm2fasta.pl";
my $CoordsReportFile = $OutPath ."/". "Coords_PresenceAbsence.csv";

#Loading presence/absence matrix
($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, @PresenceAbsenceMatrix) = Matrix($PresenceAbsence);

for ($i=1;$i<$LinesOnPresenceAbsence;$i++){
    for ($j=1;$j<$ColumnsOnPresenceAbsence;$j++){
        $ORF = $PresenceAbsenceMatrix[$i][0];
        $Sample = $PresenceAbsenceMatrix[0][$j];
        chomp($ORF,$Sample);
        $Sample =~s/\s//g;
        
        $CoordsReport -> [$i][0] = $ORF;
        $CoordsReport -> [0][$j] = $Sample;
        
        $PanOrfHmm       = $ORFsPath    ."/". $ORF ."/". $ORF . ".hmm";
        $PanOrfAln       = $ORFsPath    ."/". $ORF ."/". "*" . $ORF . ".aln.fasta";
        $SampleGenome    = $GenomesPath ."/". $Sample . ".fasta";
        $OrfOutDir       = $OutPath     ."/". "ORFs/". $ORF;
        #$TempOutDir      = $OutPath     ."/". $Temp;
        $HmmHits         = $OutPath     ."/". $Sample ."_". $ORF . ".txt";
        $SampleGeneFile  = $OrfOutDir   ."/". $Sample ."_". $ORF . ".fasta";
        $TempAln         = $OrfOutDir   ."/". $Sample ."_". $ORF . ".temp.aln.sto";
        $CoreOrfAln      = $OrfOutDir   ."/". $ORF . ".aln.sto";
        $CoreOrfHmm      = $OrfOutDir   ."/". $ORF . ".hmm";
        
        MakeDir($OrfOutDir);
        #MakeDir($TempOutDir);

        @MultifastaSampleGenome = ReadMultiFastaFile($SampleGenome);
        foreach $SingleSeq(@MultifastaSampleGenome){
            ($Header,$Seq) = ReadSeq($SingleSeq);
            $Seq{$Header} = $Seq;
        }
        
        if($j==1){
            $Hmm = $PanOrfHmm;
            $Aln = $PanOrfAln;
        }else{
            $Hmm = $CoreOrfHmm;
            $Aln = $CoreOrfAln;
        }
        
        system("nhmmer -E $eValue --noali --cpu 20 --tblout $HmmHits $Hmm $SampleGenome");
        
        @HmmHits = ReadFile($HmmHits);
        $BestHit = $HmmHits[0];
        $BestHit =~ s/\s^//g;
        $BestHit =~ s/\s+/,/g;
        @BestHitData = split(",",$BestHit);
        chomp@BestHitData;
        $ContigName = $BestHitData[0];
        $Coord1     = $BestHitData[6];
        $Coord2     = $BestHitData[7];
        
        if($Coord2 > $Coord1){
            $GeneStart = $Coord1-1;
            $GeneEnd   = $Coord2;
        }else{
            $GeneStart = $Coord2-1;
            $GeneEnd   = $Coord1;
        }
        $GeneLen = $GeneEnd-$GeneStart;
        
        $CoordsReport -> [$i][$j] = $GeneStart."_".$GeneEnd;
        
        $GeneSeq = substr$Seq{$ContigName},$GeneStart,$GeneLen;
        if($Coord1 > $Coord2){
            $GeneSeq = RevCom($GeneSeq);
        }
                
        open (FILE, ">$SampleGeneFile");
        print FILE ">$Sample\n$GeneSeq\n";
        close FILE;
        
        if($j==1){
            system("cp $Aln $TempAln");
        }else{
            system("mv $Aln $TempAln");
        }
        system("hmmalign -o $CoreOrfAln --mapali $TempAln $Hmm $SampleGeneFile");
        system("rm $TempAln");
        system("hmmbuild $CoreOrfHmm $CoreOrfAln");
        system("rm $HmmHits");
    }
}

open (FILE, ">$CoordsReportFile");
for($i=0; $i<$LinesOnPresenceAbsence+1; $i++){
    for ($j=0; $j<$ColumnsOnPresenceAbsence; $j++){
        if (defined ($CoordsReport -> [$i][$j])){
            print FILE $CoordsReport -> [$i][$j], ",";      
        }else{
            $CoordsReport -> [$i][$j] = "";
            print FILE $CoordsReport -> [$i][$j], ",";
        }
    }
    print FILE "\n";
}
close FILE;
        
for ($i=1;$i<$LinesOnPresenceAbsence;$i++){
    $ORF           = $PresenceAbsenceMatrix[$i][0];
    chomp$ORF;
    $OrfOutDir     = $OutPath     ."/". "ORFs/". $ORF;      
    $CoreOrfAln    = $OrfOutDir   ."/". $ORF . ".aln.sto";
    
    my $fastaCoreOrfAln    = $OrfOutDir   ."/". $ORF . ".aln.fasta";
    my $stoCoreOrfAln = `perl $Sto2FastaScript -g $CoreOrfAln > $fastaCoreOrfAln`;
    my @CoreOrfAln = ReadMultiFastaFile($fastaCoreOrfAln);
        
    foreach my $SingleAln(@CoreOrfAln){
        ($Header, $Seq) = ReadSeq($SingleAln);
        $AlnSeq{$Header} = $Seq;
    }
    
    for ($j=1;$j<$ColumnsOnPresenceAbsence;$j++){
        $Sample        = $PresenceAbsenceMatrix[0][$j];
        chomp$Sample;
        $Sample =~s/\s//g;
                    
        $OutSeq        = $OutPath ."/". $Sample . "_Core.aln.fasta";
        $SampleOrfAln  = $OrfOutDir   ."/". $Sample ."_". $ORF . ".aln.fasta";
        
        #%AlnSeq = Sto2Fas($CoreOrfAln);
        #
        $CatSeq{$Sample} .= $AlnSeq{$Sample};
        
        open (FILE, ">>$SampleOrfAln");
        print FILE ">$Sample" .' | '. "$ORF\n$AlnSeq{$Sample}\n";
        close FILE;

        open (FILE, ">>$OutSeq");
        if ($i==1){
            print FILE ">$Sample\n$AlnSeq{$Sample}\n";
        }else{
            print FILE "$AlnSeq{$Sample}\n";
        }
        close FILE;
    }
}

exit;