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

my ($Usage, $ProjectName, $List, $TrustedORFeome, $eVal, $PIdent, $CPUs, $Help,
    $PanGenome, $CoreGenome, $Bolean, $AnnotationPath, $MolType, $Recovery);

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
        'recovery|r'     => \$Recovery
        ) or die "USAGE:\n  $0 [--help] [--project -p prefix] [--list -l filename]
      [--trusted -c filename] [--evalue -e evalue] [--ident -i percentage]
      [--cpus -t]
\n  Use \'--help\' to print detailed descriptions of options.\n\n";

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
        system("perl $FormatORFeomes $MainPath $ProjectName $List $TrustedORFeome $AnnotationPath $MolType");
        #system("perl $FilterORFeomes $ProjectName $List $TrustedORFeome $MainPath");
        system("perl $MakeBlastDb $MainPath $ProjectName $List $TrustedORFeome $MolType");
        system("perl $InitialComparison $MainPath $ProjectName $List $TrustedORFeome $MolType $eVal $PIdent $CPUs");
        system("perl $GeneContent $MainPath $ProjectName $List $MolType $eVal $CPUs $Recovery $AnnotationPath");
}elsif($Recovery == "1"){
        system("perl $GeneContent $MainPath $ProjectName $List $MolType $eVal $CPUs $Recovery $AnnotationPath");
}

system("perl $GeneContentPlot $MainPath $ProjectName $List");

if($Functions == "1"){
        $Presence = $MainPath ."/". $ProjectName ."/". $ProjectName ."/". "_Presence_Absence.csv";
        $Core     = $MainPath ."/". $ProjectName ."/". $ProjectName ."/". "_CoreGenome.csv";
        system("perl $Functions $MainPath $ProjectName $AnnotationPath $Presence 1 1");
        system("perl $Functions $MainPath $ProjectName $AnnotationPath $Core 0 1");
}

if($PanGenome){
        system("perl $ConsensusSeq $MainPath $ProjectName 1 0");
}

if($CoreGenome){
        system("perl $CoreAlign $MainPath $ProjectName 125");
}

if($Bolean){
        system("perl $BoleanPresenceAbsence $MainPath $ProjectName");
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