#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>
#include "sislin.h"
#include "utils.h"
#include <string.h>

// Alocaçao de matriz em memória.
SistLinear_t *alocaSisLin(unsigned int n)
{
  SistLinear_t *SL = (SistLinear_t *)malloc(sizeof(SistLinear_t));

  if (SL)
  {

    SL->n = n;
    SL->A = (double **)malloc(n * sizeof(double *));
    SL->b = (double *)malloc(n * sizeof(double));

    if (!(SL->A) || !(SL->b))
    {
      liberaSisLin(SL);
      return NULL;
    }

    // Matriz como vetor de N ponteiros para um único vetor com N*N elementos
    SL->A[0] = (double *)malloc(n * n * sizeof(double));
    if (!(SL->A[0]))
    {
      liberaSisLin(SL);
      return NULL;
    }

    for (int i = 1; i < n; ++i)
    {
      SL->A[i] = SL->A[i - 1] + n;
    }
  }

  return (SL);
}

// Liberacao de memória
void liberaSisLin(SistLinear_t *SL)
{
  if (SL)
  {
    if (SL->A)
    {
      if (SL->A[0])
        free(SL->A[0]);
      free(SL->A);
    }

    if (SL->b)
      free(SL->b);

    free(SL);
  }
}

// /***********************
//  * Função que gera os coeficientes de um sistema linear k-diagonal
//  * i,j: coordenadas do elemento a ser calculado (0<=i,j<n)
//  * k: numero de diagonais da matriz A
//  ***********************/
static inline double generateRandomA(unsigned int i, unsigned int j, unsigned int k)
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return ((i == j) ? (double)(k << 1) : 1.0) * (double)rand() * invRandMax;
}

// /***********************
//  * Função que gera os termos independentes de um sistema linear k-diagonal
//  * k: numero de diagonais da matriz A
//  ***********************/
static inline double generateRandomB(unsigned int k)
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return (double)(k << 2) * (double)rand() * invRandMax;
}

/*!
  \brief Cria coeficientes e termos independentes do SL
  *
  \param SL Ponteiro para o sistema linear
  \param tipo Tipo de sistema linear a ser criado. Pode ser:
     comSolucao, eqNula, eqProporcional, eqCombLinear, hilbert
  \param coef_max Maior valor para coeficientes e termos independentes
*/
void iniSisLin(SistLinear_t *SL, unsigned int nDiagonais)
{
  for (int i = 0; i < SL->n; ++i)
  {
    for (int j = 0; j < SL->n; ++j)
    {
      if (i == j)
        SL->A[i][j] = generateRandomA(i, j, nDiagonais);
      else if (i > j)
      {

        int resp = j + (nDiagonais / 2);
        if (resp >= i)
        {
          SL->A[i][j] = generateRandomA(i, j, nDiagonais);
        }
        else
          SL->A[i][j] = 0.0;
      }
      else
      {

        int resp = j - (nDiagonais / 2);
        if (resp <= i)
        {
          SL->A[i][j] = generateRandomA(i, j, nDiagonais);
        }
        else
          SL->A[i][j] = 0.0;
      }
    }
    
    SL->b[i] = generateRandomB(nDiagonais);
  }
}

/***********************************************************************/

void prnSisLin(SistLinear_t *SL)
{
  int n = SL->n;

  for (int i = 0; i < n; ++i)
  {
    printf("\n  ");
    for (int j = 0; j < n; ++j)
      printf("%10g", SL->A[i][j]);
    printf("   |   %g", SL->b[i]);
  }
  printf("\n\n");
}

void prnVetor(double *v, unsigned int n)
{
  int i;

  printf("\n");
  for (i = 0; i < n; ++i)
    printf("%10g ", v[i]);
  printf("\n\n");
}

/**
 * @brief - Printa o vetor V no arquivo arqSaida
 * 
 * @param v - vetor
 * @param n - tamanho do vetor
 * @param arqSaida - arquivo de saida
 */
void prnVetorArq(double *v, unsigned int n, FILE *arqSaida)
{
  int i;
  fprintf(arqSaida, "\n");
  for (i = 0; i < n; ++i)
  {
    fprintf(arqSaida, "%.15g ", v[i]);
  }
  fprintf(arqSaida, "\n");
}

void prnMat(double **mat, unsigned int n, unsigned int m)
{
  for (unsigned int i = 0; i < n; ++i)
  {
    for (unsigned int j = 0; j < m; ++j)
      printf(" %f", mat[i][j]);
    printf("\n");
  }
}

/**
 * @brief - Calcula o produto entre dois vetores que geram um valor
 * produto = vetA + vetB
 * @param vetA - Vetor A
 * @param vetB - Vetor B
 * @param n - Tamanho do vetor
 * @return double - Produto dos vetores
 */
double multiplicaVetores(double *vetA, double *vetB, unsigned int n)
{
  double produto = 0;

  for (int i = 0; i < n; ++i)
  {
    produto = produto + vetA[i] * vetB[i];
    // Testa por valores inválidos.
    if (isnan(produto) || isinf(produto))
    {
      fprintf(stderr, "Erro variavel invalida: produto(multiplicaVetores): %g é NaN ou +/-Infinito\n", produto);
      exit(1);
    }
  }
  return produto;
}

/**
 * @brief Copia vetor
 *
 * @param a vetor origem
 * @param b vetor destino
 * @param n tamanho do vetor
 */
void copiaVetor(double *a, double *b, unsigned int n)
{
  int i;

  for (i = 0; i < n; i++)
  {
    b[i] = a[i];
  }
}
