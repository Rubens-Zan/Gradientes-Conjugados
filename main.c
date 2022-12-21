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
	fprintf(arqSaida,"###########\n");

    double *matIdentidade = (double *) malloc(sizeof(double)* SL->n); 

    iniSisLin(SL, comando->nDiagonais);
    precondicionador_identidade(SL,matIdentidade);
    aplicaPreCondicSL(SL, matIdentidade); 

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

    gradienteConjugadoPreCondic(SL, matIdentidade, comando->nIter,comando->erroMax,matSaida);

	fclose(arqSaida);

    free(comando);
    free(matIdentidade);
    free(x);

    liberaSisLin(SL);
    return 0;
}