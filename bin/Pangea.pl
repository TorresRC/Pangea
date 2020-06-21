#!/usr/bin/perl -w
#################################################################################
#                                                                               #
# Programador:   Roberto C. Torres, Alfonso Méndez Tenorio                      #
# Fecha:         July 2nd 2017                                                  #
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

my ($List, $eVal, $PIdent, $CPUs, $Help, $PanGenome, 
    $CoreGenome, $Boolean, $AnnotationPath, $MolType, $OutPath,
    $Recovery, $PCore);

$Recovery = 0;
$CPUs = 2;
$PCore = 99;
GetOptions(
        'help'           => \$Help,
        'list|l=s'       => \$List,
        'evalue|e=f'     => \$eVal,
        'ident|i=i'      => \$PIdent,
        'cpus|t=i'       => \$CPUs,
        'conpan|P'       => \$PanGenome,
        'alncore|C'      => \$CoreGenome,
        'pcore|c=s'      => \$PCore,
        'booleantbl|b'   => \$Boolean,
        'annotation|a=s' => \$AnnotationPath,
        'moltype|m=s'    => \$MolType,
        'out|o=s'        => \$OutPath,
        'recovery|r'     => \$Recovery
        ) or die $Usage;

if($Help){
        print "
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

print "\nList: $List\ne: $eVal\nPI: $PIdent\nCPUs $CPUs\nMol type: $MolType\n\n";

my ($FormatORFeomes, $FilterORFeomes, $MakeBlastDb, $InitialComparison,
    $PresenceAbsenceTbl, $GeneContent, $GeneContentPlot, $BooleanPresenceAbsence,
    $ConsensusSeq, $CoreAlign, $Start, $End, $Time, $RunTime, $Period,
    $ORFsFunctions, $PresenceAbsenceAnnotation, $PresenceAbsenceFile, $Core,
    $ProgressFile, $GeneContentRScript, $JobsInParallel, $FormatORFeomesStep,
    $MakeBlastDbStep, $InitialComparisonStep, $CoreFeaturesAlignStep, $Query,
    $Reference, $DataFiles, $FeaturesAlign);
my ($i);
my (@List);

$Start = time();

$FormatORFeomes            = $Src ."/". "FormatORFeomes.pl";
#$FilterORFeomes            = $Src ."/". "FilterORFeomes.pl";
$MakeBlastDb               = $Src ."/". "MakeBlastDBs.pl";
$InitialComparison         = $Src ."/". "InitialComparison.pl";
$PresenceAbsenceTbl        = $Src ."/". "PresenceAbsenceTbl.pl";
$FeaturesAlign             = $Src ."/". "FeaturesAlign.pl";

$GeneContent               = $Src ."/". "GeneContent.pl";
$GeneContentPlot           = $Src ."/". "GeneContentPlot.pl";
$BooleanPresenceAbsence    = $Src ."/". "BooleanPresenceAbsence.pl";
$ConsensusSeq              = $Src ."/". "ConsensusSeq.pl";
$CoreAlign                 = $Src ."/". "CoreAlign.pl";
$ORFsFunctions             = $Src ."/". "ORFsFunctions.pl";
$PresenceAbsenceAnnotation = $Src ."/". "PresenceAbsenceAnnotation.pl";

$ProgressFile          = $OutPath ."/". 'Progress.csv';
$GeneContentPlot       = $OutPath ."/". "GeneContentPlot.pdf";
$GeneContentRScript    = $OutPath ."/". "GeneContentScript.R";
$PresenceAbsenceFile   = $OutPath ."/". "Presence_Absence.csv";
$Core                  = $OutPath ."/". "CoreGenome.csv";

$DataFiles             = $OutPath   ."/". "DataFiles";
$FormatORFeomesStep    = $DataFiles ."/". "FormattingStep.temp";
$MakeBlastDbStep       = $DataFiles ."/". "BlastDBsBuildingStep.temp";
$InitialComparisonStep = $DataFiles ."/". "InitialComparisonStep.temp";
$CoreFeaturesAlignStep = $DataFiles ."/". "CoreFeaturesAlignStep.temp";

MakeDir($OutPath);
MakeDir($DataFiles);
@List = ReadFile($List);

$JobsInParallel = $CPUs;

if($Recovery == "0"){
    for ($i=0;$i<@List;$i++){
        system("echo perl $FormatORFeomes $List[$i] $AnnotationPath $MolType $OutPath | cat >> $FormatORFeomesStep");
        system("echo perl $MakeBlastDb $List[$i] $MolType $OutPath | cat >> $MakeBlastDbStep");
        if ($i > 0){
            system("echo perl $InitialComparison $List[$i] $List[0] $MolType $eVal $PIdent $CPUs $OutPath | cat >> $InitialComparisonStep");
        }
    }
        
    print "Formatting annotated features...";
    system("cat $FormatORFeomesStep | parallel -j$JobsInParallel");
    print "Done!";
        
    print "\nBuilding BLAST databases...";
    system("cat $MakeBlastDbStep | parallel -j$JobsInParallel");
    print "Done!";
        
    print "\nMake initial comparisons...";
    system("cat $InitialComparisonStep | parallel -j$JobsInParallel");
    print "Done!\n";
    
    system("perl $PresenceAbsenceTbl $List $PCore $OutPath");
    
    system("perl $FeaturesAlign $PCore $OutPath $CoreFeaturesAlignStep");
    print "\nAligning core features...";
    system("cat $CoreFeaturesAlignStep | parallel -j$JobsInParallel");
    print "Done!\n";
    
    
        
        #system("perl $GeneContent $ProjectName $List $MolType $eVal $CPUs $OutPath $Recovery $AnnotationPath");
}elsif($Recovery == "1"){
        #system("perl $GeneContent $ProjectName $List $MolType $eVal $CPUs $OutPath $Recovery $AnnotationPath");
}

#system("perl $GeneContentPlot $ProgressFile $ProjectName $List $OutPath");
#system("perl $ORFsFunctions $ProjectName $AnnotationPath $PresenceAbsenceFile $OutPath");
#system("perl $PresenceAbsenceAnnotation $ProjectName PanGenes $AnnotationPath $PresenceAbsenceFile $OutPath");
#system("perl $PresenceAbsenceAnnotation $ProjectName CoreGenes $AnnotationPath $Core $OutPath");
#
#if($PanGenome){
#        system("perl $ConsensusSeq $ProjectName PanGenome $PresenceAbsenceFile $OutPath/ORFs $OutPath");
#}
#
#if($CoreGenome){
#        system("perl $CoreAlign $ProjectName 125 $Core $OutPath/ORFs $OutPath");
#}
#
#if($Boolean){
#        system("perl $BooleanPresenceAbsence $ProjectName $PresenceAbsenceFile $OutPath");
#}

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

print "\n\tFinished! The analysis took $RunTime $Period\n\n";
exit;
