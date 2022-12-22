#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>
#include "utils.h"

double timestamp(void)
{
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
  return((double)(tp.tv_sec + tp.tv_nsec*1.0e-9));
}


/**
 * @brief Função para fazer o tratamento da entrada
 * @param argc 
 * @param argv 
 * @param comando 
 */
void tratamentoEntrada(int argc, char **argv, tComando *comando){
    comando->erroMax = 0; // se o erro for 0, vai ser desconsiderado
    comando->nIter = 0;
    comando->dimensao = 0;
    comando->dimensao = 0;
    
    for (int i=0;i<argc;++i){
        if (strcmp ( argv[i], "-n") == 0 && (i+1 < argc)){
            comando->dimensao = atoi(argv[i+1]); // dimensao do sl
            i++; 
        }else if (strcmp ( argv[i], "-k") == 0 && (i+1 < argc)){
           comando->nDiagonais = atoi(argv[i+1]); // numero de diagonais
           i++;
        }else if (strcmp ( argv[i], "-p") == 0 && (i+1 < argc)){
            comando->usarPreCondicionador = argv[i+1] > 0; // sem precondicionador ou de jacobi 
            i++;
        }else if(strcmp ( argv[i], "-i") == 0 && (i+1 < argc)){
            comando->nIter = atoi(argv[i+1]);  // numero de iteraçoes
            i++;
        }else if(strcmp ( argv[i], "-e") == 0 && (i+1 < argc)){
            comando->erroMax = atoi(argv[i+1]);  // erro maximo
            i++;           
        }else if (strcmp ( argv[i], "-o") == 0 && (i+1 < argc)){
           strcpy(comando->saida, argv[i+1]); //arquivo de saida
           i++;
        }  
    }

    if (!comando->dimensao || !comando->nDiagonais || !comando->nIter){
        fprintf (stderr, "ERRO argumento invalido\n");
		exit(1);
    }

}


double ** alocarMatriz(int lin,int col){
    double **matriz;
    matriz = malloc(lin * sizeof(double*)); 
    matriz[0]=malloc(lin *col*sizeof(double));

    for (int i=1; i < lin; i++)
        matriz[i] = matriz[0] + i * col;

    return matriz; 
}


void liberarMatriz(double **matriz){
    free(matriz[0]);
    free(matriz);
}

void inicializarMatriz(double **vet, int lin,int col){
    for (int i=0;i < lin;++i)
        for (int j=0;j < col;++j){
            vet[i][j]=0; 
        }

}
