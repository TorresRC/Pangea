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
#Split file by comma separated values 
sub SplitCSV{
    my ($Row) = @_;
    my @SplitedRow = split(',',$Row);
    chomp @SplitedRow;
    
    return @SplitedRow;
}

################################################################################
#Split file by tab separated values 
sub SplitTab{
    my ($Row) = @_;
    my @SplitedRow = split('\t',$Row);
    chomp @SplitedRow;
    
    return @SplitedRow;
}

################################################################################
#Split file by space separated values 
sub SplitSpaced{
    my ($Row) = @_;
    my @SplitedRow = split('\s',$Row);
    chomp @SplitedRow;
    
    return @SplitedRow;
}

################################################################################
#Split file by chart separated values 
sub SplitCharted{
    my ($Row) = @_;
    my @SplitedRow = split('\[',$Row);
    chomp @SplitedRow;
    
    return @SplitedRow;
}

################################################################################
#Esta rutina verifica si un archivo existe o no
#$name contiene la ruta completa al archivo que analizamos
sub SearchFile{
    my ($name) = @_;
    if (-e $name) {
        return 1;   #Si el archivo existe da 1   
    } else {
        return 0;   #Si no existe da 0
    }
}

################################################################################
sub SearchDir{
    my ($Dir) = @_;
    if (-d "$Dir"){
    }else{
    	my $cmd = `mkdir $Dir`;
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
#ReadFile subrutine reads a whole file and puts it in an array
sub DataFromFile{
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
			if ($Row =~/^\w/) {
				push @File, $Row;
			}	
        }
        return @File;
}

#################################################################################
sub ReadMultiFastaFile{
    my ($InputFile) = @_;
    
    $/=">";       

    open (FILE, $InputFile);      
        my $HeaderChar = <FILE>;
        my @Seq = <FILE>;
        chomp @Seq;
    close FILE;
    
    $/="\n";
    
    return @Seq;
}

#################################################################################

sub ReadSeq{
    my ($InputSeq) = @_;
    my ($Seq, @SingleFasta);
    my ($Header, @Seq) = split('\n', $InputSeq);
    chomp ($Header, @Seq);
    $Header =~ s/\n//g;
    $Header =~ s/\s//g;
    $Seq = join('',@Seq);
    $Seq =~ s/\n//g;
    $Seq =~ s/\s//g;
    $Seq =~ tr/acgt/ACGT/;
    #my @OutSeq = split('',$Seq);

    return ($Header, $Seq);
}

################################################################################
#ReadSeq subroutine reads a sequence from a fasta file and returns it in the $Seq
#variable
sub ReadSeqFile{
        my ($SeqFileName) = @_;
        my ($SeqTitle, $Seq);
        my (@Seq, @DataSeq);
        unless (open (FILE, $SeqFileName)){
                print "The Routine ReadSeq can not open $SeqFileName file\n\tExit!\n";
                exit;
        }
       $SeqTitle = <FILE>;
       chomp $SeqTitle;
       @Seq = <FILE>;
       chomp @Seq;
       close FILE;
       
       $Seq = join('',@Seq);
       $Seq =~ s/\n//g;
       $Seq =~ s/\s//g;
       $Seq =~ tr/acgt/ACGT/;
       
        return $Seq;
}

################################################################################
#ReadSeq subroutine reads a sequence from a fasta file and returns it in the $Seq
#variable
sub ReadcSeqFile{
        my ($SeqFileName) = @_;
        my ($Title, $Seq);
        my @Seq;
        unless (open (FILE, $SeqFileName)){
                print "The Routine ReadSeq can not open $SeqFileName file\n\tExit!\n";
                exit;
        }
       $Title = <FILE>;
       chomp $Title;
       @Seq = <FILE>;
       chomp @Seq;
       close FILE;
       
       $Seq = join('',@Seq);
       $Seq = reverse$Seq;
       $Seq =~ s/\n//g;
       $Seq =~ s/\s//g;
       $Seq =~ tr/acgt/ACGT/;
       $Seq =~ tr/ACGT/TGCA/;
       
        return $Seq;
}

################################################################################

sub ReadPredictFile{
    my ($InputFile) = @_;
    unless(open(FILE, $InputFile)){
        print "ReadPredict subrutine cant read $InputFile\n\tExit!";
        exit;
    }
    my @File = <FILE>;
    foreach my $Row (@File){
                if ($Row =~/^>/) {
                   shift @File;     
                }        
        }
        close FILE;       
        return @File;    
}

################################################################################

sub ReadHmmTbl{
    my ($InputFile) = @_;
    unless(open(FILE, $InputFile)){
        print "ReadHmmTbl subrutine cant read $InputFile\n\tExit!";
        exit;
    }
    my @File = <FILE>;
    my @Data;
    foreach my $Row (@File){
                if ($Row =~/^#.*/){
                    #next;
                }else{
#                    shift @File;
                    push @Data, $Row;
                }
    
    }
    close FILE;       
    return @Data;    
}


1;