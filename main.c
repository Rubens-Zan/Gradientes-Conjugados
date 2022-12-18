#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"

int main(int argc, char **argv){
    tComando *comando = (tComando *)malloc(sizeof(tComando));
    tratamentoEntrada(argc,argv,comando);
    SistLinear_t * SL = alocaSisLin(comando->dimensao);
    iniSisLin(SL, comando->nDiagonais);
    
    prnSisLin(SL); 

    // gradienteConjugado();
    free(comando); 
    liberaSisLin(SL); 
    return 0;
}