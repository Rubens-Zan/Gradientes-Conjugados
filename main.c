#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"
#include "resolvedorGradConjug.h"


int main(int argc, char **argv)
{
    srand(20222);
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

    if (comando->usarPreCondicionador){
        gradienteConjugadoPreCondic(SL, comando->nIter,comando->erroMax,arqSaida);
    } else {
        gradienteConjugado(SL,comando->nIter,comando->erroMax, arqSaida); 
    }

	fclose(arqSaida);
    free(comando);
    liberaSisLin(SL);
    return 0;
}