#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"

int main(int argc, char **argv)
{
    tComando *comando = (tComando *)malloc(sizeof(tComando));
    tratamentoEntrada(argc, argv, comando);
    SistLinear_t *SL = alocaSisLin(comando->dimensao);
    iniSisLin(SL, comando->nDiagonais);

    // prnSisLin(SL);

    // TESTANDO FUNCOES
    int n = 2;
    double **resid = alocarMatriz(2, n + 1);
    double **residAnt = alocarMatriz(2, n + 1);
    double **A = alocarMatriz(3, 3);
    double **proxX = alocarMatriz(n+1, 2);
    double **xAnt = alocarMatriz(n+1, 2);
    double **p = alocarMatriz(n+1, 2);

    
    xAnt[0][0] = 0.2356;
    xAnt[1][0] = 0.3384;

    p[0][0] = -0.3511;
    p[0][1] = -0.7229;


    resid[0][0] = -0.2810;
    resid[0][1] = 0.7492;
    residAnt[0][0] = -8;
    residAnt[0][1] = -3;
    A[0][0] = 4;
    A[0][1] = 1;
    A[1][0] = 1;
    A[1][1] = 3;
    // printf("beta: %f \n", calcBeta(resid, residAnt, n)); OK
    // printf("alpha: %f \n", calcAlpha(residAnt, A, residAnt, n)); OK
    calcProxX(proxX,xAnt,0.4122, p, n);
    prnMat(proxX, n,1); 

    //////

    // gradienteConjugado();
    free(comando);
    liberaSisLin(SL);
    return 0;
}