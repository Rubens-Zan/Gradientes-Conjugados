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
    
    //TESTANDO FUNCOES 
    int n =2; 
    double resid[2][n+1] = {{-0.2810, 0.7492}};
    double residAnt[2][n+1] = {{-8, -3}};
    double A[2+1][2+1] = {{4,1},{1,3}}; 

    printf("beta: %f \n",calcBeta(resid,residAnt,n));
    printf("alpha: %f \n",calcAlpha(residAnt,A, residAnt, n));
    
    //////

    // gradienteConjugado();
    free(comando); 
    liberaSisLin(SL); 
    return 0;
}