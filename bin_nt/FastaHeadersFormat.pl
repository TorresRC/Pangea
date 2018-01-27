#/usr/bin/perl -w
#This script changes each sequence's header in a file. The header icludes a counter
#in a xxn format
#By Roberto Torres
use strict;
use FindBin;
use lib "$FindBin::Bin/../lib_nt";
use Routines;

my ($MainPath, $InputFile, $SeqPath, $i, $n, $FileName, $Prefix, $OutFile, $Row, $Count,
	$Counter, $OutHeader, $Ext, $OutPath);
my (@List, @Assembly);

$MainPath = "/Users/bioinformatica/Documents/CoreGenome";

print "\nPlease type the name of the file list which incloudes all sequences you want to edit: ";
$InputFile = <STDIN>;
chomp $InputFile;
#$InputFile = $MainPath ."/". "Metadata/GenomesUsedForCore.ls";

print "\nPlease type the full path where the query sequences are located: ";
$SeqPath = <STDIN>;
chomp $SeqPath;
#$SeqPath = $MainPath ."/". "Sequences/genomes";

#print "\nPlease type the file extension that you want to save the outputs: ";
#$Ext = <STDIN>;
#chomp $Ext;
$Ext = ".ffn";

$OutPath = $SeqPath ."/". "HeadersFormated";
MakeDir($OutPath);

@List = ReadFile($InputFile);

$n = scalar@List;
for ($i=0;$i<$n;$i++){
	$FileName = $List[$i];
	$Prefix = Prefix($FileName);
	
	$InputFile = $SeqPath ."/". $FileName .".ffn";
	@Assembly = ReadFile($InputFile);

	$OutFile = $OutPath."/".$Prefix.$Ext;
	
	open(FILE, ">>$OutFile");
	$Count = 0;
	foreach $Row (@Assembly){	
		if ($Row =~ "^>" ){
			$Counter = Counter($Count);
			my @RowData = split (" ", $Row);
			$OutHeader = $RowData[0];
			chomp $OutHeader;
			#$OutHeader = ">".$Prefix."_Contig".$Counter;
			#$OutHeader = ">Contig".$Counter;
			#$OutHeader = ">$Prefix" . "_$Counter";
			#$OutHeader = ">".$Prefix;
			#$OutHeader = ">Gene" . "_$Counter";
			print FILE "$OutHeader\n";
			$Count = $Counter;
		}else{
			print FILE "$Row\n";
		}
	}
	Progress($n, $i);
}
exit;
