#!/usr/bin/perl -w
#################################################################################
#                                                                               #
# Programador:   Roberto C. Torres, Alfonso Méndez Tenorio                      #
# Fecha:         2 de julio de 2017                                             #
#################################################################################
use strict;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my $Src      = "$FindBin::Bin";
my $MainPath = "$FindBin::Bin/../..";
my $Usage    = "USAGE:\n  $0 [--help] [--project -p prefix] [--list -l filename]
      [--trusted -c filename] [--evalue -e evalue] [--ident -i percentage]
      [--cpus -t] [--out -o path]
\n  Use \'--help\' to print detailed descriptions of options.\n\n";

my ($ProjectName, $List, $TrustedORFeome, $eVal, $PIdent, $CPUs, $Help,
    $PanGenome, $CoreGenome, $Bolean, $AnnotationPath, $MolType, $OutPath,
    $Recovery);

$Recovery = 0;
GetOptions(
        'help'           => \$Help,
        'project|p=s'    => \$ProjectName,
        'list|l=s'       => \$List,
        'trusted|c=s'    => \$TrustedORFeome,
        'evalue|e=f'     => \$eVal,
        'ident|i=i'      => \$PIdent,
        'cpus|t=i'       => \$CPUs,
        'conpan|P'       => \$PanGenome,
        'alncore|C'      => \$CoreGenome,
        'boleantbl|b'    => \$Bolean,
        'annotation|a=s' => \$AnnotationPath,
        'moltype|m=s'    => \$MolType,
        'out|o=s'        => \$OutPath,
        'recovery|r'     => \$Recovery
        ) or die $Usage;

if($Help){
        print "
        \t--project      <Project_Name>
        \t--list         <List_File_Name>
        \t--trusted      <TrustedFile.fasta>
        \t--evalue       <e-value>
        \t--ident        <Percetage_of_Identity>
        \t--cpus         <CPUs>
        \t--conpan       <Bolean>
        \t--alncore      <Bolean>
        \t--boleantbl    <Bolean>
        \t--annotatedtbl <Bolean>
        \t--moltype      <nucl or prot>
        \t--recovery     <Bolean>
        \n\n";
        exit;
}

print "\nProject: $ProjectName\nList: $List\nTrusted: $TrustedORFeome\ne: $eVal\nPI: $PIdent\nCPUs $CPUs\nMol type: $MolType\n\n";

my ($Project, $FormatORFeomes, $SortGenes, $FilterORFeomes, $MakeBlastDb,
    $InitialComparison, $GeneContent, $GeneContentPlot, $BoleanPresenceAbsence,
    $ConsensusSeq, $CoreAlign, $Start, $End, $Time, $RunTime, $Period,
    $Functions, $Presence, $Core);

$Start = time();

$FormatORFeomes        = $Src ."/". "FormatORFeomes.pl";
$SortGenes             = $Src ."/". "SortGenes.pl";
$FilterORFeomes        = $Src ."/". "FilterORFeomes.pl";
$MakeBlastDb           = $Src ."/". "MakeBlastDBs.pl";
$InitialComparison     = $Src ."/". "InitialComparison.pl";
$GeneContent           = $Src ."/". "GeneContent.pl";
$GeneContentPlot       = $Src ."/". "GeneContentPlot.pl";
$BoleanPresenceAbsence = $Src ."/". "BoleanPresenceAbsence.pl";
$ConsensusSeq          = $Src ."/". "ConsensusSeq.pl";
$CoreAlign             = $Src ."/". "CoreAlign.pl";
$Functions             = $Src ."/". "ORFsFunctions.pl";

if($Recovery == "0"){
        system("perl $FormatORFeomes $ProjectName $List $TrustedORFeome $AnnotationPath $MolType $OutPath");
        #system("perl $FilterORFeomes $ProjectName $List $TrustedORFeome $MainPath");
        system("perl $MakeBlastDb $ProjectName $List $TrustedORFeome $MolType $OutPath");
        system("perl $InitialComparison $ProjectName $List $TrustedORFeome $MolType $eVal $PIdent $CPUs $OutPath");
        system("perl $GeneContent $ProjectName $List $MolType $eVal $CPUs $OutPath $Recovery $AnnotationPath");
}elsif($Recovery == "1"){
        system("perl $GeneContent $ProjectName $List $MolType $eVal $CPUs $OutPath $Recovery $AnnotationPath");
}

system("perl $GeneContentPlot $ProjectName $List $OutPath");

$Presence = $OutPath ."/". $ProjectName ."/". $ProjectName ."/". "_Presence_Absence.csv";
$Core     = $OutPath ."/". $ProjectName ."/". $ProjectName ."/". "_CoreGenome.csv";
system("perl $Functions $ProjectName $AnnotationPath $Presence $OutPath 1 1");
system("perl $Functions $ProjectName $AnnotationPath $Core $OutPath 0 1");

if($PanGenome){
        system("perl $ConsensusSeq $ProjectName $OutPath 1 0");
}

if($CoreGenome){
        system("perl $CoreAlign $ProjectName $OutPath 125");
}

if($Bolean){
        system("perl $BoleanPresenceAbsence $ProjectName $OutPath");
}

#Timestamp
$End = time();
$Time = ($End-$Start);
if ($Time < 3600){
        $RunTime = ($Time)/60;
        $Period = "minutes";
}elsif ($Time >= 3600 && $Time < 86400){
        $RunTime = (($Time)/60)/60;
        $Period = "hours";
}elsif ($Time >= 86400){
        $RunTime = ((($Time)/60)/60)/24;
        $Period = "days";
}

print "\n\tFinished! The $ProjectName gene content analysis took $RunTime $Period\n\n";
exit;