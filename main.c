#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"
#include "opmatrizes.h"


int main(int argc, char **argv)
{
    tComando *comando = (tComando *)malloc(sizeof(tComando));
    tratamentoEntrada(argc, argv, comando);
    // SistLinear_t *SL = alocaSisLin(comando->dimensao+1);
    SistLinear_t *SL = alocaSisLin(2);
    FILE *arqSaida;
	double matSaida[comando->nIter+1][2];

	arqSaida = fopen(comando->saida,"w+");
	fprintf(arqSaida,"###########\n");
    double **x = alocarMatriz(SL->n+1,2); 
    double **matPreConj = alocarMatriz(SL->n+1,SL->n+1); 
    geraMatrizIdentidade(matPreConj,SL->n); 

    // iniSisLin(SL, comando->nDiagonais);
    SL->A[0][0] = 4;
    SL->A[0][1] = 1;
    SL->A[1][0] = 1;
    SL->A[1][1] = 3;

    SL->b[0] = 1;
    SL->b[1] = 2;
    x[0][0] = 2;
    x[1][0] = 1;

    // inicializarMatriz(x, SL->n, 1); // inicializa chute inicial com 0 
    prnSisLin(SL); 
    printf("\n"); 
    gradienteConjugado(SL,x,matPreConj,comando->nIter,comando->erroMax, matSaida);


	fclose(arqSaida);

    free(comando);
    liberaSisLin(SL);
    liberarMatriz(x); 
    liberarMatriz(matPreConj); 
    return 0;
}