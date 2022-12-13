#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>
#include "funcAux.h"

// double timestamp(void)
// {
//   struct timespec tp;
//   clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
//   return((double)(tp.tv_sec + tp.tv_nsec*1.0e-9));
// }

// /***********************
//  * Função que gera os coeficientes de um sistema linear k-diagonal
//  * i,j: coordenadas do elemento a ser calculado (0<=i,j<n)
//  * k: numero de diagonais da matriz A
//  ***********************/
inline double generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return ( (i==j)?(double)(k<<1) : 1.0 )  * (double)rand() * invRandMax;
}

// /***********************
//  * Função que gera os termos independentes de um sistema linear k-diagonal
//  * k: numero de diagonais da matriz A
//  ***********************/
inline double generateRandomB( unsigned int k )
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return (double)(k<<2) * (double)rand() * invRandMax;
}

/**
 * @brief Função para fazer o tratamento da entrada
 * @param argc 
 * @param argv 
 * @param comando 
 */
void tratamentoEntrada(int argc, char **argv, tComando *comando){
    for (int i=0;i<argc;++i){
        if (strcmp ( argv[i], "-n") == 0 && (i+1 < argc)){
            comando->dimensao = atoi(argv[i+1]); 
            i++; 
        }else if (strcmp ( argv[i], "-k") == 0 && (i+1 < argc)){
           comando->nDiagonais = atoi(argv[i+1]);
           i++;
        }else if (strcmp ( argv[i], "-p") == 0 && (i+1 < argc)){
            comando->usarPreCondicionador= argv[i+1] > 0; 
            i++;
        }else if(strcmp ( argv[i], "-i") == 0 && (i+1 < argc)){
            comando->nIter = atoi(argv[i+1]); 
            i++;
        }else if(strcmp ( argv[i], "-e") == 0 && (i+1 < argc)){
            comando->erroMax = atoi(argv[i+1]); 
            i++;           
        }else if (strcmp ( argv[i], "-o") == 0 && (i+1 < argc)){
           strcpy(comando->saida, argv[i+1]);
           i++;
        }  
    }

  printf("dimensao:%d err:%d ndiago:%d n iter:%d precond:%d saida: %s\n",comando->dimensao,comando->erroMax,comando->nDiagonais,comando->nIter,comando->usarPreCondicionador, comando->saida); 
}