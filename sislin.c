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

/***********************************************************************/

/**
 * @brief Calcula a proxima direcao de busca p<i>
 * p<i> = r<i> + beta<i-1> * p<i-1> 
 * @param resid - Residuo calculado r<i>
 * @param beta - Beta calculado beta<i-1>
 * @param direcAnterior - Direcao anterior calculada 
 * @return double 
 */
double **calcProxDirecBusca(double **resid, double beta,double **direcAnterior, int n){
    double proxDir[n+1][n+1];
    for (int i=0; i < n;++i)
        for(int j=0;j < n;++j)
            proxDir[i][j] = resid[i][j] + beta * direcAnterior[i][j]; 
    return proxDir; 
}


/**
 * @brief Calcula beta com base no residuo atual e no anterior
 * 
 * @param resid 
 * @param residAnt 
 * @return double 
 */
double calcBeta(double *resid,double *residAnt, int n){
    
}

/**
 * @brief Calcula o proximo x<i>
 * x<i> = x<i-1> + a <i> * p<i-1> 
 * @param xAnt - Valor do x<i-1> anterior
 * @param a - valor do a<i-1> anterior
 * @param p - valor da direção de busca p<i-1)> anterior
 * @return double - Proximo x calculado
 */
double **calcProxX(double **xAnt, double a, double **p, int n){
    double xAtual[n+1][n+1];
    for (int i =0;i < n;++i)
        for (int j=0;j<n;++j)
            xAtual[i][j] = xAnt[i][j] + a * p[i][j]; 
    return xAtual; 
}


/**
 * @brief Metodo resolvedor de sistemas lineares pelo metodo de gradiente conjugado
 * 
 * @param A  - MATRIZ A SER TESTADA
 * @param b - Coeficiente de solucao
 * @param x -  Solucao do sistema
 * @param M  - Matriz precondicionadora 
 * @param maxIt - NUMERO MAXIMO DE ITERACOES
 * @param tol - TOLERANCIA
 * @param n  - DIMENSAO 
 * @return int - NUMERO DE ITERACOES
 */
int gradienteConjugado(double *A, double *b, double *x, double *M, int maxIt, double tol, int n){

}

/***********************************************************************/

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


