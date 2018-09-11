#!/usr/bin/perl -w
#################################################################################
#Scipt ORFsFunction.pl                                                          #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          10 de enero de 2018                                             #
#################################################################################
use strict; 
use List::MoreUtils qw{any first_index};
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $ProjectName, $AnnotationPath, $OutPath, $PresenceAbsence,
    $ORFsFunctionsList, $AnnotatedTable);

$Usage = "\tUsage: ORFsFunctions.pl <Main_Path> <Project Name> <Annotation_Path> <Presence-Absence>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName       = $ARGV[0];
$AnnotationPath    = $ARGV[1];
$PresenceAbsence   = $ARGV[2];
$OutPath           = $ARGV[3];
$ORFsFunctionsList = $ARGV[4];
$AnnotatedTable    = $ARGV[5];

my ($Project, $AnnotatedPresenceAnsence, $LinesOnPresenceAbsence,
    $ColumnsOnPresenceAbsence, $ORF, $cmd, $OTU, $AnnotationFile, $Function,
    $OTUORF, $Index, $ColumnsOnAnnotation, $Gene, $ECNumber, $ORFsFunctionsFile,
    $Prefix, $Line, $LocusTag
    );
my ($i,$j);
my (@PresenceAbsenceMatrix, @Annotation, @PresenceAbsence, @Header, @ORFData,
    @InFileName, @AnnotationFile);
my (%Annotation, %RepresentantLocus, %RepresentantOTU, %LocusTag, %Function, %RepresentantFunction, %Gene, %ECN,
    %RepresentantGene, %RepresentantECN);
my $Annotated = [ ];

@InFileName = split("/",$PresenceAbsence);
$Prefix = $InFileName[3];

$AnnotatedPresenceAnsence = $OutPath ."/". $ProjectName . $Prefix . "_Annotated_Presence_Absence.csv";
$ORFsFunctionsFile        = $OutPath ."/". $ProjectName . $Prefix . "_ORFsAnnotation.csv";
        
open (FILE, ">$ORFsFunctionsFile");
    print FILE "ORF,Gene,EC_Number,Product";
close FILE;

print "\nLoading Presence Absence File...";
($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, @PresenceAbsenceMatrix) = Matrix($PresenceAbsence);
print "Done!";
        
print "\nFinding Representant OTUs:\n";
for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
    $ORF = $PresenceAbsenceMatrix[$i]->[0];
    for ($j=1; $j< $ColumnsOnPresenceAbsence; $j++){
        $OTU = $PresenceAbsenceMatrix[0]->[$j];
        $OTU =~s/\r//g;
        $LocusTag{$ORF}{$OTU} = $PresenceAbsenceMatrix[$i]->[$j];
        if ($LocusTag{$ORF}{$OTU} ne ""){
            $RepresentantLocus{$ORF} = $LocusTag{$ORF}{$OTU};
            $RepresentantOTU{$ORF} = $OTU;
            last;
        }else{
            next;
        }
    }
    Progress($LinesOnPresenceAbsence,$i);
}

print "\nObtaining Functions:\n";  
for ($i=1; $i<$ColumnsOnPresenceAbsence; $i++){
    $OTU = $PresenceAbsenceMatrix[0]->[$i];
    $OTU =~s/\r//g;
    $AnnotationFile = $AnnotationPath ."/". $OTU ."/". $OTU . ".tsv";
    @AnnotationFile = ReadFile($AnnotationFile);
             
    foreach $Line (@AnnotationFile){
        @Annotation = split("\t",$Line);
        $LocusTag = $Annotation[0];
        if(scalar@Annotation == 3){
            $Gene = "";
            $ECNumber = "";
        }elsif(scalar@Annotation == 4){
            if ($Annotation[2] =~ /^\d/){
                $ECNumber = $Annotation[2];
                chomp$ECNumber;
                $Gene = "";
            }else{
                $Gene = $Annotation[2];
                chomp$Gene;
                $ECNumber = "";
            }
        }elsif(scalar@Annotation == 5){
            $Gene = $Annotation[2];
            $ECNumber = $Annotation[3];
        }      
        $Function = $Annotation[$#Annotation];
        chomp$Function;
        $Function =~ s/,/-/g;
        $Gene{$LocusTag}{$OTU}     = $Gene;
        $ECN{$LocusTag}{$OTU}      = $ECNumber;
        $Function{$LocusTag}{$OTU} = $Function;  
    }
    Progress($ColumnsOnPresenceAbsence,$i);
}

print "\nAnnotating ORFs:\n";
for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
    $ORF = $PresenceAbsenceMatrix[$i]->[0];
        $RepresentantGene{$ORF}     = $Gene{$RepresentantLocus{$ORF}}{$RepresentantOTU{$ORF}};
        $RepresentantECN{$ORF}      = $ECN{$RepresentantLocus{$ORF}}{$RepresentantOTU{$ORF}};
        $RepresentantFunction{$ORF} = $Function{$RepresentantLocus{$ORF}}{$RepresentantOTU{$ORF}};
    Progress($LinesOnPresenceAbsence,$i);
}

print "\nWriting File:\n";
for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
    $ORF = $PresenceAbsenceMatrix[$i]->[0];
        open (FILE, ">>$ORFsFunctionsFile");
            print FILE "\n$ORF,$RepresentantGene{$ORF},$RepresentantECN{$ORF},$RepresentantFunction{$ORF}";
        close FILE;
    Progress($LinesOnPresenceAbsence,$i);
}
exit;        

#if ($AnnotatedTable == "1"){
#        ($LinesOnPresenceAbsence, $ColumnsOnPresenceAbsence, @PresenceAbsenceMatrix) = Matrix($PresenceAbsence);
#        
#        $Annotated -> [0][0] = $PresenceAbsenceMatrix[0]->[0];
#        
#        print "\nGetting annotation of each ORF:\n";
#        for ($i=1; $i<$LinesOnPresenceAbsence; $i++){
#           for ($j=1; $j<$ColumnsOnPresenceAbsence; $j++){
#              $OTU = $PresenceAbsenceMatrix[0][$j];
#              $OTU =~s/\r//g;
#        
#              $Annotated -> [$i][0] = $PresenceAbsenceMatrix[$i]->[0];
#              $Annotated -> [0][$j] = $OTU;
#              $ORF = $PresenceAbsenceMatrix[$i][$j];
#              
#              $AnnotationFile = $AnnotationPath ."/". $OTU ."/". $OTU . ".tsv";
#             
#              if ($ORF ne ""){
#                 $cmd = `grep \"$ORF\tCDS\" $AnnotationFile`;
#                 
#                 @Annotation = split("\t",$cmd);
#                 $Function = $Annotation[$#Annotation];
#                 $Function =~ s/,/-/g;
#                 $Function =~ s/\b//g;
#                 $Function =~ s/ /_/g;
#                 $Function =~ s/\n//g;
#                 $Function =~ s/\s/_/g;
#                 $Function =~ s/\//-/g;
#                 chomp$Function;
#                 $Annotated -> [$i][$j] = $Function;
#        
#              }else{
#                 $Annotated -> [$i][$j] = "";
#              }
#           }
#           Progress($LinesOnPresenceAbsence, $i);
#        }
#        
#        print "\nBuilding Annotated file:\n";
#        open (FILE, ">$AnnotatedPresenceAnsence");
#        for ($i=0; $i<$LinesOnPresenceAbsence; $i++){
#           for ($j=0; $j<$ColumnsOnPresenceAbsence; $j++){
#              print FILE $Annotated -> [$i][$j], ",";
#           }
#           print FILE "\n";
#           Progress($LinesOnPresenceAbsence, $i);
#        }
#        close FILE;
#}
exit;