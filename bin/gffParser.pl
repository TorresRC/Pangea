   #!/usr/bin/perl -w
#################################################################################
#By:            Roberto C. Torres, PhD                                          #
#e-mail:        torres.roberto.c@gmail.com                                      #
#Date:          June 10th, 2020                                                 #
#################################################################################

$GffFile $ARGV[0];


$/='##';
unless(open (FILE, $GffFile)){
  print "error\nCan not open $GffFile file called on $0\n";
  exit;
}
@GffFile = <FILE>;
chomp @GffFile;
close FILE;
$/="\n";
chomp @GffFile;

foreach $Block (@GffFile){
  @Block = split('\n',$Block);
  $BlockHeader = shift@Block;
  if(scalar@Block > 1){
    if($BlockHeader eq 'FASTA'){                        #gff fasta block
      @Contigs = split('>',$Block);
      shift@Contigs;
      foreach $Contig (@Contigs){
        ($ContigName,@ContigSeq) = split("\n",$Contig);
        push @{$ContigsOfSample{$Sample}},$ContigName;
        $ContigSeq{$Sample}{$ContigName}  = join("",@ContigSeq);
        $ContigLength{$Sample}{$ContigName} = length$ContigSeq{$Sample}{$ContigName};
      }
    }else{                                         #gff annotation block
      foreach $Line (@Block){
      @GeneInfo   = split("\t",$Line);
      @GeneContext = split (';',$GeneInfo[8]);
      $Gene       = $GeneContext[0];
      $Gene       =~ s/^ID=//g;
      $Function   = $GeneContext[$#GeneContext];
      $Function   =~ s/^product=//g;
      $ContigName = $GeneInfo[0];
      $Contig{$Sample}{$Gene} = $ContigName;
      $StartCoord{$Sample}{$ContigName}{$Gene} = $GeneInfo[3];
      $EndCoord{$Sample}{$ContigName}{$Gene}   = $GeneInfo[4];
      $Strand{$Sample}{$ContigName}{$Gene}     = $GeneInfo[6];
      $Function{$Sample}{$ContigName}{$Gene}   = $Function;
    }
  }  
}
