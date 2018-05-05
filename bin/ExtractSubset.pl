#!/usr/bin/perl -w
#################################################################################
#Scipt ConsensusPanGenome.pl                                                    #
#                                                                               #
#Programmer:    Roberto C. Torres                                               #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          19 de octubre de 2017                                           #
#################################################################################
use strict; 
use List::MoreUtils qw{any first_index};
use FindBin;
use lib "$FindBin::Bin/../lib";
use Routines;

my ($Usage, $PresenceAbsence, $List, $SubsetFile);

$Usage = "\tUsage: ConsensusSeq.pl <Original Presence/Absence File> <Subset Presence/Absence>\n";
unless(@ARGV) {
        print $Usage;
        exit;
}
chomp @ARGV;
$PresenceAbsence = $ARGV[0];
$List  = $ARGV[1];
$SubsetFile = $ARGV[2];

my($ORF, $OTU, $Qry, $LinesOnPresenceabsenceFile, $ColumnsOnPresenceabsenceFile);
my($i, $j, $k, $nSubset, $nColumns, $SubsetColumns, $Index, $Gene, $Strain);
my(@List, @PresenceAbsenceMatrix);
my %GeneOfStrain;
my $Subset = [ ];

print "Loading the Presence/Absence file...";
($LinesOnPresenceabsenceFile, $ColumnsOnPresenceabsenceFile, @PresenceAbsenceMatrix) = Matrix($PresenceAbsence);
print "Done!\n";

for ($i=1;$i<$LinesOnPresenceabsenceFile;$i++){
  for ($j=1;$j<$ColumnsOnPresenceabsenceFile;$j++){
   $ORF = $PresenceAbsenceMatrix[$i]->[0];
   $Strain = $PresenceAbsenceMatrix[0]->[$j];
   $Strain =~ s/\s//g;
   $Strain =~ s/\n//g;
   $ORF =~ s/\s//g;
   $ORF =~ s/\n//g;
   
   if ($PresenceAbsenceMatrix[$i]->[$j] ne ""){
      $Gene = $PresenceAbsenceMatrix[$i]->[$j];
   }else{
      $Gene = "";
   }
   $GeneOfStrain{$ORF}{$Strain} = $Gene;
  }
}

@List = ReadFile($List);
$nSubset = scalar@List;

$Subset -> [0][0] = "ORF";

print "Building a Subset presence/absence file:\n";
for ($i=1;$i<$LinesOnPresenceabsenceFile;$i++){
   for ($j=0;$j<$nSubset;$j++){
      $ORF     = $PresenceAbsenceMatrix[$i]->[0];
      $Strain  = $List[$j];
      $Strain =~ s/\s//g;
      $ORF =~ s/\s//g;
      $Subset -> [$i][0] = $ORF;
      $Subset -> [0][$j+1] = $Strain;
      $Subset -> [$i][$j+1] = $GeneOfStrain{$ORF}{$Strain};
   }
}
  
open (FILE, ">$SubsetFile");
for ($i=0; $i<$LinesOnPresenceabsenceFile; $i++){
     for ($j=0; $j<$nSubset+1; $j++){
      $Gene = $Subset -> [$i][$j];
        if ($j < $nSubset){
                print FILE $Gene, ",";
        }else{
                print FILE $Gene;
        }
     }
     print FILE "\n";
     Progress($LinesOnPresenceabsenceFile, $i);
}
close FILE;

exit;