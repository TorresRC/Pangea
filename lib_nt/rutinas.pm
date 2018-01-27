#################################################################################
# Modulo:        Coleccion de rutinas para la manipulacion de secuencias de     #
#                DNA.                                                           #
# Programador:   Alfonso Méndez Tenorio                                         #
#                                                                               #
# Fecha:         24 de febrero de 2016                                          #
# Actualización:                                                                #
#     14/III/2017: Roberto incluyó subrutina MakeDir                            #
#################################################################################
use strict;

#Esta subrrutina muestra el contenido de archivos de un directorio
#En este caso el directorio contiene archivos con extension .fasta
#$Path contiene la ruta completa del directorio que analizamos
sub MuestraDir {
        my ($Path) = @_;
        opendir(my $dh, $Path) || die;
        while(readdir $dh) {
            if (($_ ne '.') && ($_ ne '..') && ($_=~ /\.fasta/)){
                print "$_\n";
            }
        }
        closedir $dh;
        #Esta rutina no devuelve valores solo imprime los archivos
}

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

# Esta subrutina lee una secuencia de DNA desde el archivo FASTA de nombre
# $Nombrearchivo.
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
# mas conveniente para almacenar esta informacion pero es conveniente
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

#Esta rutina busca el directorio especificado, si no existe éste es creado 
sub MakeDir{
    my ($NewDir) = @_;
	if (-d "$NewDir"){
	}else{
		my $cmd = `mkdir $NewDir`;
	}
}

#ReadFile subrutine reads a file and puts it in an array
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
                }else{ push @File, $Row;     
                }        
        }
        return @File;
}

#Prefix subrutine set a prefix from the name of each file
sub Prefix{
        my ($FileName) = @_;
        my @SplitName = split ('\.',$FileName);
        my $Prefix = $SplitName[0];
        my $Ext = $SplitName[1];
        
        return $Prefix;    
}

1;
