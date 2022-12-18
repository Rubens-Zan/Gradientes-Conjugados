#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"
#include "opmatrizes.h"


int main(int argc, char **argv)
{
    tComando *comando = (tComando *)malloc(sizeof(tComando));
    tratamentoEntrada(argc, argv, comando);
    SistLinear_t *SL = alocaSisLin(comando->dimensao+1);
    double **x = alocarMatriz(SL->n+1,2); 
    double **matPreConj = alocarMatriz(SL->n+1,SL->n+1); 
    geraMatrizIdentidade(matPreConj,SL->n); 

    iniSisLin(SL, comando->nDiagonais);
    
    gradienteConjugado(SL,x,matPreConj,comando->nIter,comando->erroMax);

    free(comando);
    liberaSisLin(SL);
    liberarMatriz(x); 
    liberarMatriz(matPreConj); 
    return 0;
}