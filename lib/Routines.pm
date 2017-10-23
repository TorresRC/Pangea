#!/usr/bin/perl -w
################################################################################
# Routine colection to DNA parsing                                             #
# By Roberto Torres                                                            #
################################################################################
use strict;

#################################################################################

#################################################################################
sub Progress{
    my ($n, $i, $null) = @_;
    $i++;
    my $Percentage = ($i/$n)*100;
    
    if ($i<$n){
      my $Progress = sprintf "%.3d", $Percentage;
      print "\r\tProgress: [$Progress%]";
      $|=1;
    }else{
      print "\r\tProgress: [100%]\n\n";
      $|=1;
    }
}

#################################################################################
sub Counter{
    my ($Count) = @_;
    $Count++;
    my $Counter = sprintf "%.4d", $Count;
    return $Counter;  
}

################################################################################
sub MakeDir{
    my ($NewDir) = @_;
	if (-d "$NewDir"){
	}else{
		my $cmd = `mkdir $NewDir`;
        return $cmd;
	}
}

################################################################################
#Prefix subrutine set a prefix from the name of each file
sub Prefix{
        my ($FileName) = @_;
        my @SplitName = split ('\.',$FileName);
        my $Prefix = $SplitName[0];
        my $Ext = $SplitName[1];
        
        return $Prefix;    
}

################################################################################
#ReadFile subrutine reads a whole file and puts it in an array
sub ReadFile{
        my ($InputFile) = @_;
        unless (open(FILE, $InputFile)){
               print "The Routine ReadFile Can not open $InputFile file\n";
               exit;
        }
        my @Temp = <FILE>;
        chomp @Temp;
        close FILE;
        my @File;
        foreach my $Row (@Temp){
               if ($Row =~/^#/) {
               }else{
				push @File, $Row;     
               }
        }
        return @File;
}
################################################################################
sub AnnotatedGenes{
        my ($File) = @_;
        my $cmd = `grep ">" $File`;
           $cmd =~ s/>//g;
           $cmd =~ s/\h//g;
        my @Data = split('\n',$cmd);
        return @Data;
}
################################################################################

sub GenesInBlastReport{
        my ($File, $GeneId, $null) = @_;
        open (FILE, ">>$File");
                print FILE "$GeneId\n";
        close FILE;
}
################################################################################

sub DismissORFs{
        my ($Id, @IDs, $null) = @_;
        my $n = scalar@IDs;
        for(my $i=0;$i<$n;$i++){
                if($IDs[$i] eq $Id){
                        splice @IDs, $i, 1;
                        $n--;
                }
        }
        return @IDs;
}
################################################################################

sub Extract{
        my ($Qry, $DataBase,$Entry,$OutSeq, $null) = @_;
        print "\tExtracting ORF from $Qry...";	
        my $cmd = `blastdbcmd -db $DataBase -dbtype nucl -entry "$Entry" -out $OutSeq`;
        print "Done!\n";
}
################################################################################

sub Align{
        my ($Seq1, $Seq2, $ToAlign, $AlnFile, $null) = @_;
        print "\tAligning sequences...";
        my $cmd = `cat $Seq1 $Seq2 > $ToAlign`;
        $cmd = `muscle -in $ToAlign -out $AlnFile -quiet`;
        print "Done!\n";
}
################################################################################

sub HMM{
        my ($CPUs, $HmmFile, $AlnFile, $null) = @_;
        print "\tBuilding a HMM...";
        my $cmd = `hmmbuild --dna --cpu $CPUs $HmmFile $AlnFile`;
        print "Done!\n";
}

1;