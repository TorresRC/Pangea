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
      [--evalue -e evalue] [--ident -i percentage] [--cpus -t] [--out -o path]
\n  Use \'--help\' to print detailed descriptions of options.\n\n";

my ($ProjectName, $List, $eVal, $PIdent, $CPUs, $Help, $PanGenome, 
    $CoreGenome, $Boolean, $AnnotationPath, $MolType, $OutPath,
    $Recovery);

$Recovery = 0;
GetOptions(
        'help'           => \$Help,
        'project|p=s'    => \$ProjectName,
        'list|l=s'       => \$List,
        'evalue|e=f'     => \$eVal,
        'ident|i=i'      => \$PIdent,
        'cpus|t=i'       => \$CPUs,
        'conpan|P'       => \$PanGenome,
        'alncore|C'      => \$CoreGenome,
        'booleantbl|b'   => \$Boolean,
        'annotation|a=s' => \$AnnotationPath,
        'moltype|m=s'    => \$MolType,
        'out|o=s'        => \$OutPath,
        'recovery|r'     => \$Recovery
        ) or die $Usage;

if($Help){
        print "
        \t--project      <Project_Name>
        \t--list         <List_File_Name>
        \t--evalue       <e-value>
        \t--ident        <Percetage_of_Identity>
        \t--cpus         <CPUs>
        \t--conpan       <Boolean>
        \t--alncore      <Boolean>
        \t--booleantbl   <Boolean>
        \t--annotatedtbl <Boolean>
        \t--moltype      <nucl or prot>
        \t--recovery     <Boolean>
        \n\n";
        exit;
}

print "\nProject: $ProjectName\nList: $List\nTrusted: $TrustedORFeome\ne: $eVal\nPI: $PIdent\nCPUs $CPUs\nMol type: $MolType\n\n";

my ($Project, $FormatORFeomes, $SortGenes, $FilterORFeomes, $MakeBlastDb,
    $InitialComparison, $GeneContent, $GeneContentPlot, $BooleanPresenceAbsence,
    $ConsensusSeq, $CoreAlign, $Start, $End, $Time, $RunTime, $Period,
    $ORFsFunctions, $PresenceAbsenceAnnotation, $PresenceAbsenceFile, $Core,
    $ProgressFile, $GeneContentRScript);

$Start = time();

$FormatORFeomes            = $Src ."/". "FormatORFeomes.pl";
$SortGenes                 = $Src ."/". "SortGenes.pl";
$FilterORFeomes            = $Src ."/". "FilterORFeomes.pl";
$MakeBlastDb               = $Src ."/". "MakeBlastDBs.pl";
$InitialComparison         = $Src ."/". "InitialComparison.pl";
$GeneContent               = $Src ."/". "GeneContent.pl";
$GeneContentPlot           = $Src ."/". "GeneContentPlot.pl";
$BooleanPresenceAbsence    = $Src ."/". "BooleanPresenceAbsence.pl";
$ConsensusSeq              = $Src ."/". "ConsensusSeq.pl";
$CoreAlign                 = $Src ."/". "CoreAlign.pl";
$ORFsFunctions             = $Src ."/". "ORFsFunctions.pl";
$PresenceAbsenceAnnotation = $Src ."/". "PresenceAbsenceAnnotation.pl";

$ProgressFile          = $OutPath ."/". $ProjectName . '_Progress.csv';
$GeneContentPlot       = $OutPath ."/". $ProjectName . "_GeneContentPlot.pdf";
$GeneContentRScript    = $OutPath ."/". $ProjectName . "_GeneContentScript.R";
$PresenceAbsenceFile   = $OutPath ."/". $ProjectName . "_Presence_Absence.csv";
$Core                  = $OutPath ."/". $ProjectName . "_CoreGenome.csv";

if($Recovery == "0"){
        system("perl $FormatORFeomes $ProjectName $List $AnnotationPath $MolType $OutPath");
        system("perl $MakeBlastDb $ProjectName $List $MolType $OutPath");
        system("perl $InitialComparison $ProjectName $List $MolType $eVal $PIdent $CPUs $OutPath");
        system("perl $GeneContent $ProjectName $List $MolType $eVal $CPUs $OutPath $Recovery $AnnotationPath");
}elsif($Recovery == "1"){
        system("perl $GeneContent $ProjectName $List $MolType $eVal $CPUs $OutPath $Recovery $AnnotationPath");
}

system("perl $GeneContentPlot $ProgressFile $ProjectName $List $OutPath");
system("perl $ORFsFunctions $ProjectName $AnnotationPath $PresenceAbsenceFile $OutPath");
system("perl $PresenceAbsenceAnnotation $ProjectName PanGenes $AnnotationPath $PresenceAbsenceFile $OutPath");
system("perl $PresenceAbsenceAnnotation $ProjectName CoreGenes $AnnotationPath $Core $OutPath");

if($PanGenome){
        system("perl $ConsensusSeq $ProjectName PanGenome $PresenceAbsenceFile $OutPath/ORFs $OutPath");
}

if($CoreGenome){
        system("perl $CoreAlign $ProjectName 125 $Core $OutPath/ORFs $OutPath");
}

if($Boolean){
        system("perl $BooleanPresenceAbsence $ProjectName $PresenceAbsenceFile $OutPath");
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
