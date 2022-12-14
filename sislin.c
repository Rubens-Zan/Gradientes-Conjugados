#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "sislin.h"

// Alocaçao de matriz em memória. 
SistLinear_t* alocaSisLin (unsigned int n)
{
  SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));
  
  if ( SL ) {
    
    SL->n = n;
    SL->A = (real_t **) malloc(n * sizeof(real_t *));
    SL->b = (real_t *) malloc(n * sizeof(real_t));

    if (!(SL->A) || !(SL->b)) {
      liberaSisLin(SL);
      return NULL;
    }

    // Matriz como vetor de N ponteiros para um único vetor com N*N elementos
    SL->A[0] = (real_t *) malloc(n * n * sizeof(real_t));
    if (!(SL->A[0])) {
      liberaSisLin(SL);
      return NULL;
    }
    
    for (int i=1; i < n; ++i) {
      SL->A[i] = SL->A[i-1]+n;
    }
  }
  
  return (SL);
}

// Liberacao de memória
void liberaSisLin (SistLinear_t *SL)
{
  if (SL) {
    if (SL->A) {
      if (SL->A[0]) free (SL->A[0]);
    free(SL->A);
    }
    
    if (SL->b) free(SL->b);

    free(SL);
  }
}


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

/*!
  \brief Cria coeficientes e termos independentes do SL
  *
  \param SL Ponteiro para o sistema linear
  \param tipo Tipo de sistema linear a ser criado. Pode ser: 
     comSolucao, eqNula, eqProporcional, eqCombLinear, hilbert 
  \param coef_max Maior valor para coeficientes e termos independentes
*/
void iniSisLin (SistLinear_t *SL, unsigned int nDiagonais){
  unsigned int n = SL->n;
  // inicializando sequencia de numeros aleatorios
  srand(20222);

  // inicializa vetor b
  for (unsigned int i=0; i<n; ++i) {
    SL->b[i] = generateRandomB(nDiagonais);
  }
    
  // inicializa a matriz A
  for (unsigned int i=0; i<n; ++i) {
    for (unsigned int j=0; j<n; ++j)  {
      SL->A[i][j] = generateRandomA(i,j,nDiagonais);
    }
  }
  
}


void prnSisLin (SistLinear_t *SL)
{
  int n=SL->n;

  for(int i=0; i < n; ++i) {
    printf("\n  ");
    for(int j=0; j < n; ++j)
      printf ("%10g", SL->A[i][j]);
    printf ("   |   %g", SL->b[i]);
  }
  printf("\n\n");
}

void prnVetor (real_t *v, unsigned int n)
{
  int i;

  printf ("\n");
  for(i=0; i < n; ++i)
      printf ("%10g ", v[i]);
  printf ("\n\n");

}

