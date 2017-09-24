#!/usr/bin/perl -w

#################################################################################
#   Programa Filter Features                                                    #
#   Nota: Se debe ajustar la ruta de lib y de la variable $PathSeq a las que    #
#   realmente se tengan donde se instalción el programa.                        #
#                                                                               #
# Programador:   Roberto C. Torres                                              #
# Fecha:         10 de enero de 2017                                            #
#################################################################################

use strict; 
use lib '/home/bioinformatica/CoreGenome/src/lib';
#use lib '/Users/Roberto/Documents/lib';
use rutinas;


my ($Usage, $ProjectName, $List);

$Usage = "\tUsage: FilterFeatures.pl <Project Name> <List File Name>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$ProjectName = $ARGV[0];
$List = $ARGV[1];

my ($MainPath, $Project, $FeaturesPath, $InPath, $OutPath, $ExcludedPath,
    $MainList, $ext, $Qry, $ptt, $File, $Prefix, $OutFile, $Header, $Row,
    $Excluded, $LogFile);
my (@List, @Body);

#$MainPath = "/Users/Roberto/CoreGenome";
$MainPath = "/home/bioinformatica/CoreGenome";
$Project = $MainPath ."/". $ProjectName;
#$FeaturesPath = $Project."/". "Features";
$FeaturesPath = $MainPath ."/". "Features";
$InPath = $FeaturesPath ."/". "ptt";
$OutPath = $InPath ."/". 'Filtered';
$ExcludedPath = $OutPath ."/". "ExcludedFeatures";
$MainList = $Project ."/". $List;
$ext = ".ptt";
$LogFile = $Project ."/". $ProjectName . ".log";


open (STDERR, "| tee -ai $LogFile") or die "$0: dup: $!";

MakeDir($OutPath);
MakeDir($ExcludedPath);

@List = ReadFile($MainList);
chomp @List;

foreach $ptt (@List){
    print "\nFiltering $ptt...";
    $File = $ptt . $ext;
    $Qry = $InPath ."/". $File;
	$OutFile = $OutPath ."/Filtered_". $File;
    $Excluded = $ExcludedPath ."/". "Excluded_" . $File;
    
    unless (open(FILE, $Qry)){
                print "Can't open $Qry file\n\tFinished\n";
                exit;
    }
        $Header = <FILE>;
        @Body = <FILE>;
        chomp $Header;
        chomp @Body;
    close FILE;
    
    open (FILE, ">$OutFile");
        print FILE "$Header\n";
    close FILE;
    
    foreach $Row (@Body){
        if ($Row =~ /^>/){
            open (FILE, ">>$Excluded");
            print FILE "$Row\n";
            close FILE;
        }else{
            open (FILE, ">>$OutFile");
            print FILE "$Row\n";
            close FILE;
        }
    }
    print "Done.";
}
print "\n\n\tFinished!\n\n";
exit;
