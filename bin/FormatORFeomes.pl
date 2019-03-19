#################################################################################
#Scipt FormatORFeomes.pl                                                        #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          19 de octubre de 2017                                           #
#################################################################################
use strict;
use FindBin;
use lib "$FindBin::Bin/../lib";
use List::MoreUtils qw{any};
use Routines;

my ($Usage, $ProjectName, $List, $TrustedName, $AnnotationPath, $MolType, $OutPath);

$Usage = "\tUsage: FormatORFeomes.pl <Main_Path> <Project Name> <List File Name> <Annotation Path> <Molecule Type>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;

$ProjectName    = $ARGV[0];
$List           = $ARGV[1];
$TrustedName    = $ARGV[2];
$AnnotationPath = $ARGV[3];
$MolType        = $ARGV[4];    #nucl or prot
$OutPath        = $ARGV[5];

my ($Project, $MainList, $ORFeomesPath, $Ext, $Qry, $Prefix, $OriginalORFeome,
	 $FormatedORFeome, $Line, $FormatedHeader);
my ($i, $n);
my (@List, @ORFeome, @ORFHeader);

#$Project      = $OutPath ."/". $ProjectName;
$Project      = $OutPath;
$ORFeomesPath = $Project ."/". "ORFeomes";

if ($MolType eq "nucl"){
	$Ext = ".ffn";
}elsif($MolType eq "prot"){
	$Ext = ".faa";
}

MakeDir($OutPath);
MakeDir($ORFeomesPath);

@List = ReadFile($List);
if ( any {$_ ne $TrustedName} @List){
        unshift @List, $TrustedName;
}
$n = scalar@List;

for ($i=0;$i<$n;$i++){
	$Qry = $List[$i];
	$Prefix = Prefix($Qry);
	$OriginalORFeome = $AnnotationPath ."/". $Prefix ."/". $Prefix.$Ext;
	
	@ORFeome = ReadFile($OriginalORFeome);
	$FormatedORFeome = $ORFeomesPath ."/". $Prefix.$Ext;	
	open(FILE, ">>$FormatedORFeome");
	foreach $Line (@ORFeome){	
		if ($Line =~ "^>" ){
			@ORFHeader = split (" ", $Line);
			$FormatedHeader = $ORFHeader[0];
			chomp $FormatedHeader;
			print FILE "$FormatedHeader\n";
		}else{
			$Line =~ s/-//g;
			print FILE "$Line\n";
		}
	}
	Progress($n, $i);
}
exit;