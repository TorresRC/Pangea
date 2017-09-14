################################################################################
# Modulo:        Colección de rutinas para la manipulación de secuencias de    #
#                DNA.                                                          #
# Programador:   Alfonso Méndez Tenorio                                        #
# Fecha:         24 de agosto de 2004                                          #
# Actualización:                                                               #
# 9-Marzo-2006  Se adicionaron rutinas para la búsqueda, fuerza bruta y el     #
#               de Knut-Morris-Pratt (KMP).                                    #
#               Los algoritmos fueron traducidos del Pascal a Perl en la forma #
#               como se describen en el libro "Introduction to Algorithms in   #
#               Pascal" de T. W. Parsons.                                      #
#               Hay diferencias importantes con el código en Pascal ya que en  #
#               en este lenguaje los caracteres de los strings se numeran   -  #
#               desde 1, pero en Perl se numeran desde 0.                      #
#               Se incluye una versión del algoritmo KMP pero sin la tabla     #
#               esto hace que el KMP se comporte entonces como el brute search.#
#               La tabla del KMP completo es una FSA y esta es la clave de la  #
#               eficiencia de este algoritmo que nunca retrocede y por lo tanto#
#               realiza la búsqueda con menos pasos.                           #
#               Se incluye otra versión de la búsqueda pero con expresiones    #
#               regulares. Las expresiones regulares se codifican también como #
#               FSAs pero mucho mas complejas y robustas que la del KMP.       #
################################################################################
use strict;

# Esta subrutina lee una secuencia de DNA desde el archivo FASTA de nombre
# $Nombrearchivo.
# No es la mejor rutina todavía
sub LeerDNA{
    my($Nombrearchivosub) = @_;
    my $Titulosub;
    my $DNAsub;
    my @DNAarray;
    unless (open (ARCH, $Nombrearchivosub))
    {
        print "El archivo $Nombrearchivosub no pudo abrirse\n";
        print "Programa terminado\n";
        exit;
    }
    $Titulosub = <ARCH>;
    chomp $Titulosub;
    @DNAarray = <ARCH>;
     close ARCH;
    $DNAsub = join('',@DNAarray);
    $DNAsub =~ s/\n//g;
    $DNAsub =~ s/\s//g;
    $DNAsub =~ tr/acgt/ACGT/;
    print "Titulo: $Titulosub\n";
    return $DNAsub;
}

# Esta rutina devuelve el contenido de bases de una secuencia de DNA.
# El contenido se almacena en un arreglo. Esta estructura no es la
# más conveniente para almacenar esta información pero es conveniente
# para ilustrar el manejo de los arreglos.
sub CalcularComposicion{
    my ($DNA) = @_;
    my @Composicionsub;
    my $posicion;
    my $Longitud;
    my $Base;
    $Longitud = length $DNA;
    $Composicionsub[0] = 0;
    $Composicionsub[1] = 0;
    $Composicionsub[2] = 0;
    $Composicionsub[3] = 0;
    for ($posicion=0; $posicion<$Longitud; ++$posicion) {
        $Base = substr($DNA, $posicion, 1);
        if ($Base eq 'A') {++$Composicionsub[0]};
        if ($Base eq 'C') {++$Composicionsub[1]};
        if ($Base eq 'G') {++$Composicionsub[2]};
        if ($Base eq 'T') {++$Composicionsub[3]};
    }
    return @Composicionsub;
}


# Esta subrutina calcula el valor del Tm de una molécula de DNA
# Para lo cual tomo en cuenta su composición y la longitud
sub CalcularTm {
    my ($Xgccont, $Long) = @_;
    my $Tmresult;
    $Tmresult = (81.5 + (41*($Xgccont))-(500/$Long));
    return  $Tmresult;
}

# Esta subrutina calcula el valor de la entropía de Shannon a partir
# del contenido de bases del DNA
# Se trata de una implementación muy rudimentaria. (poco eficiente)
sub CalcularEntropia {
    my($As, $Cs, $Gs, $Ts, $Long) = @_;
    my ($Htemp) = 0.0;
    #my($indice) = 0;
    my($prob)  = 0.0;
    my $indice;
    for ($indice = 0; $indice < 4; ++$indice) 
    {
        if ($indice == 0) {
            $prob = $As/$Long;  #prob es la probabilidad 
            $Htemp = $Htemp + (($prob)* log ($prob));
        }
        if ($indice == 1) {
            $prob = $Cs/$Long;
            $Htemp = $Htemp + (($prob)* log ($prob));
        }
        if ($indice == 2) {
            $prob = $Gs/$Long;
            $Htemp = $Htemp + (($prob)* log ($prob));
        }
        if ($indice == 3) {
            $prob = $Ts/$Long;
            $Htemp = $Htemp + (($prob)* log ($prob));
        }
    }
    return -$Htemp;
}


# En esta nueva versión de la rutina CalcularEntropia se acepta un arreglo con 
# el contenido de bases del DNA
# Ilustra como pasar un arreglo por referencia a una subrrutina
# Se trata de una implementación mas eficiente
sub CalcularEntropia2 {
    my($Composition) = @_;
    my ($Htemp) = 0.0;
    my($prob);
    my $indice;
    my $Long;
    $Long = @$Composition[0] + @$Composition[1] + @$Composition[2] + @$Composition[3];
    for ($indice = 0; $indice < 4; ++$indice) 
    {
        $prob = @$Composition[$indice]/$Long;
        if ($prob > 0){$Htemp = $Htemp + (($prob)* log ($prob));}
    }
    return -$Htemp;
}


#Devuelve una posición al azar en la secuencia de DNA entre 0 y Len-1
sub PosicionAleatoria{
    my ($DNAsub) = @_;
    my $Lensub;
    $Lensub = length $DNAsub;
    return int(rand($Lensub));
}


#Esta rutina pretende ilustrar el uso de las funciones para generar
#numeros aleatorios, produciendo mutaciones al azar en una secuencia de
#DNA.
#La rutina requiere del uso de otras rutinas auxiliares para aislar el proceso
#Se puede perfeccionar para que el nucleótido escogido siempre sea diferente
sub MutarDNA{
    my($DNA, $posicion) = @_;
    my @nucleotidos;
    my $Base;
    @nucleotidos = ('A', 'C', 'G', 'T');
    my $mutacion;
    $Base = substr($DNA, $posicion, 1);  #opcional
    #Seleccionar un nucleotido al azar
    do {
      $mutacion = $nucleotidos[int(rand @nucleotidos)];
    } until not ($Base eq $mutacion);  
    #print "El nucleotido selccionado es $mutacion \n";
    substr($DNA, $posicion, 1, $mutacion);
    return $DNA;
}

#Calcula el grado de identidad (similitud) entre dos secuencias de DNA de la
#misma longitud
sub CalcularSimilitud {
    my ($DNA, $DNAmutante) = @_;
    my $indice;
    my $Len;
    my $matches = 0;
    my $BaseDNA;
    my $Basemut;
    $Len = length $DNA;
    for ($indice = 0; $indice < $Len; ++$indice) {
        $BaseDNA = substr($DNA, $indice, 1);
        $Basemut = substr($DNAmutante, $indice, 1);
        if ($BaseDNA eq $Basemut) {++$matches};
    }
    return ($matches*100/$Len);
}

sub PrintMask {
    my ($DNA, $DNAmut) = @_;
    my $Len;
    my $Index;
    my $Base;
    my $Mut;
    my $Mask;
    my $TempMask = '';
    $Len = length $DNA;
    for ($Index = 0; $Index < $Len; ++$Index){
        $Base = substr($DNA, $Index, 1);
        $Mut  = substr($DNAmut, $Index, 1);
        if ($Base eq $Mut) {$Mask = ':'} else {$Mask = ' '};
        $TempMask = $TempMask . $Mask;
    }
    return $TempMask;
}

sub Traducir {
    my($codon) = @_;
    my $aa;
    my %Codigo;
    
    %Codigo = (
               'TTT' => 'F', 'TTC' => 'F', 'TTA' => 'L', 'TTG' => 'L',
               'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',
               'TAT' => 'Y', 'TAC' => 'Y', 'TAA' => '*', 'TAG' => '*',
               'TGT' => 'C', 'TGC' => 'C', 'TGA' => '*', 'TGG' => 'W',
               'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
               'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
               'CAT' => 'H', 'CAC' => 'H', 'CAA' => 'Q', 'CAG' => 'Q',
               'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
               'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I', 'ATG' => 'M',
               'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
               'AAT' => 'N', 'AAC' => 'N', 'AAA' => 'K', 'AAG' => 'K',
               'AGT' => 'S', 'AGC' => 'S', 'AGA' => 'R', 'AGG' => 'R',
               'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V',
               'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
               'GAT' => 'D', 'GAC' => 'D', 'GAA' => 'E', 'GAG' => 'E',
               'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G',
              );
    # $aa = $Codigo{$codon};
    if (exists $Codigo{$codon}){$aa = $Codigo{$codon}}
    else {$aa = '?'};
    
    return $aa;
}

sub Ribosoma {
    my ($DNA, $inicio, $fin, $Readingframe) = @_;
    my $longitud;
    my $i;
    my $codon;
    my $aa;
    my $Proteina;
    
    $Readingframe = $Readingframe-1;
    # $fin = $longitud;
    $inicio = $inicio + $Readingframe;
    $Proteina = '';
    for ($i = $inicio;$i<$fin;$i=$i+3){
         if ($fin - $i>2){
             $codon= substr($DNA,$i,3);
             $aa = Traducir($codon); 
             # print "$codon=>$aa\n";}
             $Proteina = $Proteina . $aa;
        }
    }
    return $Proteina;       
}

sub Darlimites {
    my ($longitud) = @_;
    my %limites = ();
    my $inicio;
    my $fin;
    
    do{  
        print "Indique el comienzo de la traduccion: ";
        $inicio = <STDIN>;
        chomp $inicio;
        print "Indique el final de la traduccion: ";
        $fin = <STDIN>;
        chomp $fin;
        if ($inicio>$fin){print "el final debe ser mayor que el inicio\n"}
        if ($fin>$longitud){print "el final no debe ser mayor que $longitud\n"}
    }   until ($inicio<$fin) and ($fin<$longitud); 
        $limites {'inicio'}= $inicio;
        $limites {'fin'}= $fin;
    return %limites;
}    


sub BuscaEnzima{
    
      my ($X,$Archivoenzimas) =@_;
      my @campos;
      my @Datosenzimas;
      my $i;
      my $enzima;
      my $secuencia;
      
        unless (open(AENZIMAS,$Archivoenzimas))
        {
            print"El archivo $Archivoenzimas no fue encontrado\n ";
            return ' ';      
        
        }
        @Datosenzimas = <AENZIMAS>;
        close AENZIMAS;
        foreach(@Datosenzimas)
        {
          (/REBASE/)and next;
          (/^\s*$/)and next;
          (/=/)and next;
          (/Copy/)and next;
          (/Roberts/)and next;
            @campos = split(" ",$_);
            $enzima = shift @campos;
            $secuencia = pop @campos;
            if ($enzima eq $X) {return $secuencia;}
             
        }
        return '';
   }
sub seqtoRE
{
    my ($secuencia) = @_;
    my $expre = '';
    my $longseq;
    my $i;
    my %exregcodigo=(
                    A=>'A',
                    C=>'C',
                    G=>'G',
                    T=>'T',
                    R=>'[AG]',
                    Y=>'[CT]',
                    M=>'[AC]',
                    K=>'[GT]',
                    S=>'[CG]',
                    W=>'[AT]',
                    B=>'[CGT]',
                    D=>'[AGT]',
                    H=>'[ACT]',
                    V=>'[ACG]',
                    N=>'[ACGT]');
    $secuencia =~ s/\^//g;
    $longseq = length($secuencia);
    
    for($i=0;$i<$longseq;$i++)
    {
     $expre = $expre . $exregcodigo{substr($secuencia,$i,1)}   
    }    
 return $expre;  
}

#Esta subrrutina implementa el algoritmo llamado de la fuerza bruta

sub brute_search {
    my ($textstr, $pattern, $n, $m) = @_;
    
    #En perl estas variables comienzan en 0 porque los strings se numeran desde 0
    my $pat = 0;
    my $off = 0;
    
    while (($pat < $m) and (($pat + $off) < $n)) {
        if (substr($pattern, $pat, 1) eq substr($textstr, $pat+$off, 1)){
              $pat ++
            } else {
              $off ++;
              $pat = 0;
            }
    }   
    if ($pat > ($m-1)) {return $off}
     else {return -1}
}

# Esta es una version simplificada del algoritmo KMP pero sin la tabla next
# En este caso el algoritmo funciona como el brute search
sub kmp_sin_tabla {
    my ($textstr, $pattern, $n, $m) = @_;

    my $i = 0;
    my $j = 0;
    #my @next = make_table($pattern, $m);
    while (($j < $m) and ($i<$n)){
        if ($j eq -1) {
            $i++;
            $j++;
        }else {
              if (substr($pattern, $j, 1) eq substr($textstr, $i, 1)) {
                  $i++;
                  $j++;
              } else {
                   $i = $i-$j; #retrocedemos para posicionar el patron al inicio
                               # El algoritmo kmp completo nunca retrocede.
                   $j = -1;
              } #2o. else
        } #1er else
    } #while   
    if ($j eq $m) {return $i - $m}
       else {return -1};
}

    
sub make_table{
    my ($pattern, $m) = @_;
    my @next;
    my $i = 0;
    my $j = -1;
    $next[0]=-1;
    while ($i < ($m-1)) {
        if (($j eq -1) or (substr($pattern, $i, 1) eq substr($pattern, $j, 1))) {
               $i++;
               $j++;
               if ((substr($pattern, $i, 1) ne substr($pattern, $j, 1))) {
                  $next[$i] = $j;
                } else {
                  $next[$i] = $next[$j];
                } #del primer if-else
            } else {
               $j = $next[$j];
        } #del segundo if-else
    } #del while
    return @next;
}


sub kmp_search {
    my ($textstr, $pattern, $n, $m) = @_;

    my $i = 0;
    my $j = 0;
    my @next = make_table($pattern, $m);
    while (($j < $m) and ($i<$n)){
        if ($j eq -1) {
            $i++;
            $j++;
        }else {
              if (substr($pattern, $j, 1) eq substr($textstr, $i, 1)) {
                  $i++;
                  $j++;
              } else {
                   #$i = $i-$j; 
                   $j = $next[$j];
              } #2o. else
        } #1er else
    } #while   
    if ($j eq $m) {return $i - $m}
       else {return -1};
}



1;