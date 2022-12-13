#include <stdio.h>
#include <stdlib.h>

#include "funcAux.h"

int main(int argc, char **argv){
    tComando *comando = (tComando *)malloc(sizeof(tComando));
    tratamentoEntrada(argc,argv,comando);

    free(comando); 
    return 0;
}