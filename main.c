#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"
#include "opmatrizes.h"


int main(int argc, char **argv)
{
    tComando *comando = (tComando *)malloc(sizeof(tComando));
    tratamentoEntrada(argc, argv, comando);
    SistLinear_t *SL = alocaSisLin(comando->dimensao);
    // SistLinear_t *SL = alocaSisLin(2);
    FILE *arqSaida;
	double matSaida[comando->nIter+1][2];

	arqSaida = fopen(comando->saida,"w+");
    
    
	fprintf(arqSaida,"#rzl20 Rubens Zandomenighi Laszlo \n");
	fprintf(arqSaida,"# \n");


    iniSisLin(SL, comando->nDiagonais);

    // double *matPreConj = (double *) malloc(sizeof(double)* SL->n); 
    if (comando->usarPreCondicionador){
        gradienteConjugadoPreCondic(SL, comando->nIter,comando->erroMax,matSaida,arqSaida);

    } else {
        gradienteConjugadoPreCondic(SL,comando->nIter,comando->erroMax,matSaida, arqSaida); 
    }

    // SL->A[0][0] = 4;
    // SL->A[0][1] = 1;
    // SL->A[1][0] = 1;
    // SL->A[1][1] = 3;

    // SL->b[0] = 1;
    // SL->b[1] = 2;
    // x[0][0] = 2;
    // x[1][0] = 1;

    prnSisLin(SL); 
    printf("\n"); 


	fclose(arqSaida);
    free(comando);
    liberaSisLin(SL);
    return 0;
}