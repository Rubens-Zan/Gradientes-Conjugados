#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "opmatrizes.h"
#include "sislin.h"
#include "utils.h"

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
static inline double generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return ( (i==j)?(double)(k<<1) : 1.0 )  * (double)rand() * invRandMax;
}

// /***********************
//  * Função que gera os termos independentes de um sistema linear k-diagonal
//  * k: numero de diagonais da matriz A
//  ***********************/
static inline double generateRandomB( unsigned int k )
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
// FUNCOES PARA A RESOLUCAO POR GRADIENTE CONJUGADO

/**
 * @brief Calcula a proxima direcao de busca p<k>
 * 
 * p[i][j]<k> = z[i][j]<k> + beta<k-1> * p[i][j]<k-1> ok;; 
 * 
 * @param proxDir
 * @param z - z calculado 
 * @param beta - Beta calculado beta<i-1>
 * @param direcAnterior - Direcao anterior calculada 
 * @return double 
 */
void calcProxDirecBusca(double **proxDir,double **z, double beta,double **direcAnterior, int n){
    for (int i=0; i < n;++i)
      for(int j=0;j < 1;++j)
        proxDir[i][j] = z[i][j] + beta * direcAnterior[i][j]; 
}

/**
 * @brief - Calcula alpha
 * 
 * alpha = r<k>^t*z<k> / (p<k>^t*A*p<k>)  
 * 
 * @param resid - residuo calculado
 * @param A - Matriz original
 * @param p - Matriz de direção de busca calculada
 * @param z - Matriz z
 * @param n - dimensão da matriz
 * @return double 
 */
double calcAlpha(double **resid,double **A, double **p,double **z, int n){
  double alpha = 0; 
  double **residTransp=alocarMatriz(1+1,n+1); // matriz de residuo transposta
  double **pTransp = alocarMatriz(1+1,n+1); // matriz de direcao transposta
  double **multPtxA =alocarMatriz(1+1,n+1);
  
  transporMat(resid, residTransp, n, 1);
  transporMat(p, pTransp, n, 1);
  
  printf("RESID TRANSP\n\n");
  prnMat(residTransp,1,n);

  printf("Z \n\n");

  prnMat(z,n,1);
  printf("\n\n");

  double resTxZ = formMultMatrizesGeraValor(residTransp,z,1,n,n,1); // resTxZ = resid^T * z
  printf("REST * Z \n\n", resTxZ);
  
  // multMat(pTransp,A,1,n,n,n,multPtxA); // p^T * A
  // double resPxAxP =  formMultMatrizesGeraValor(multPtxA,pTransp,1,n,n,1); // (p^T * A) * p

  // alpha = resTxZ /  resPxAxP; 

  liberarMatriz(residTransp);
  liberarMatriz(pTransp);
  liberarMatriz(multPtxA);

  return alpha; 
}

/**
 * @brief Calcula beta com base no residuo atual e no anterior
 *
 * beta<k> = (res<k+1>^T * z<k+1>)/(res<k>^T * z<k>)
 *  
 * @param resid - Residuo da k-esima iteração
 * @param residAnt - Resíduo da k-1-esima iteração
 * @return double - valor de beta calculado
 */
double calcBeta(double **resid,double **residAnt,double **z, int n){
	double beta = 0; 
  double **residTransp=alocarMatriz(1,n+1); // matriz de residuo transposta
	double **residAntTransp=alocarMatriz(1,n+1); // matriz do residuo anterior transposta
  
  transporMat(resid, residTransp, n,1);
  transporMat(residAnt, residAntTransp, n,1);

  double residTxZ = formMultMatrizesGeraValor(residTransp,z,1,n,n,1); 
  double resTantxZ = formMultMatrizesGeraValor(residAntTransp,z,1,n,n,1); 

  beta = residTxZ/ resTantxZ;

  liberarMatriz(residAntTransp);
  liberarMatriz(residTransp);

  return beta; 
}


/**
 * @brief Calcula o proximo x<k>
 *
 *  x[i][j]<k> = x[i][j]<k-1> + alpha <k> * p[i][j]<k-1>
 * 
 * @param proxX  
 * @param xAnt - Valor do x<i-1> anterior
 * @param alpha - valor do alpha<i-1> anterior
 * @param p - valor da direção de busca p<i-1)> anterior
 * @param n - dimensao do sistema linear
 * @return double - Proximo x calculado
 */
void calcX(double **proxX,double **xAnt,double alpha, double **p, int n){
    for (int i =0;i < n;++i){
      for (int j=0;j < 1;++j){
        printf("\n\n aq alpha %d \n\n");
        // proxX[i][j] = xAnt[i][j] + alpha * p[i][j]; 
      }
    }
}

/**
 * @brief - Calcula o residuo 
 * r<k> = r<k-1> - alpha<k-1> * A * p<k-1> ok
 * @param residuoAnterior 
 * @param alpha 
 * @param A 
 * @param p 
 * @param residuo 
 * @param n - dimensao do sistema linear
 */
void calcResiduo(double **residuoAnterior, double alpha, double **A, double **p, double ** residuo,int n){
  double **matAlphaxA = alocarMatriz(n+1,n+1);
  double **multAlphaxAxp = alocarMatriz(n+1,1+1); 

  for (int i=0;i < n;++i){
    for (int j=0;j<n;++j){
      matAlphaxA[i][j] = alpha * A[i][j]; // alpha * A
    }
  }

  multMat(matAlphaxA, p,n,n,n,1, multAlphaxAxp); // (alpha * A) * p  

  for (int i=0;i < n;++i)
    for(int j=0;j < 1;++j)
      residuo[i][j] = residuoAnterior[i][j]-multAlphaxAxp[i][j];  

  liberarMatriz(matAlphaxA); 
  liberarMatriz(multAlphaxAxp); 
}


/**
 * @brief Calcula o residuo inicial 
 * r0 = b - A * x0
 * @param A 
 * @param b 
 * @param x 
 * @param resid 
 */
void calcResiduoInicial(double **A,double *b,double **x,double **resid, int n){
  double **matAux = alocarMatriz(n,1); 
  multMat(A, x, n,n,n,1,matAux);
  
  for (int i=0;i < n;++i){
    for (int j=0;j < 1; ++j){
      resid[i][j] =b[i] - matAux[i][j];
    }
  }

  liberarMatriz(matAux); 
}


/**
 * @brief Metodo resolvedor de sistemas lineares pelo metodo de gradiente conjugado
 * 
 * @param SL - sistema linear
 * @param x -  Solucao do sistema inicial
 * @param M  - Matriz precondicionadora 
 * @param maxIt - NUMERO MAXIMO DE ITERACOES
 * @param tol - TOLERANCIA
 * @param matSaida - Matriz que guardara o Erro e a Norma de cada iteração
 * @return int - NUMERO DE ITERACOES
 */
int gradienteConjugadoPreCondic(SistLinear_t *SL, double **x, double **matPreConj, int maxIt, double tol, double matSaida[][2]){
    // inicia chute inicial com vetor de 0  
    double alpha, beta; 
    double **resid = alocarMatriz(SL->n+1,2); // matriz de residuo 
    double **residAnt = alocarMatriz(SL->n+1,2); // matriz de residuo anterior
    double **direc = alocarMatriz(SL->n+1,2); // matriz de direcao de busca
    double **dAnt = alocarMatriz(SL->n+1,2); // matriz de direcao anterior    
    double **xAnt = alocarMatriz(SL->n+1,2); // matriz de chute anterior
    double **z = alocarMatriz(SL->n+1,2); // matriz de 
  
    int it;
    // as inicializacoes estao ok!    
    // x = 0 
    inicializarMatriz(x,SL->n,1);  

    // r<0> = b - A * x<0>
    calcResiduoInicial(SL->A,SL->b,x,resid,SL->n);

    // z<0> = MATPRECONJ^-1 * r<0>
    multMat(matPreConj,resid,SL->n,SL->n,SL->n,1,z); 

    // inicia direcao inicial com o z inicial
    copiaMat(direc,z,SL->n, 1); // d<0> = z<0>
    
    // verifica erro do x com a tolerancia 
    // loop 
    for(it =0;it < maxIt;++it){
      // calcula alpha
      alpha = calcAlpha(resid,SL->A,direc,z,SL->n);
      printf("ALPHA: %f \n", alpha);

      copiaMat(x,xAnt,SL->n,1); // xAnterior = xAtual

      // calcula x
      // calcX(x,xAnt,alpha, direc, SL->n);
      // printf("x: \n");
      // prnMat(x,SL->n,1);

      // // calcula residuo
      // copiaMat(resid,residAnt,SL->n, 1); 
      // calcResiduo(residAnt,alpha,SL->A,direc,resid,SL->n); 
      // printf("residuo: \n");
      // prnMat(resid,SL->n,1);
      
      // // calcula z
      // // z<k+1> = C^-1 * r<k+1>
      // multMat(matPreConj,resid,SL->n,SL->n,SL->n,1,z); 
      // printf("z: \n");
      // prnMat(z,SL->n,1);
      
      // // calcula erro 
      // // ERRO = r<k+1> * r<k+1>

      // // faz o if com o  erro e tolerancia 
      
      // // calcula beta
      // beta = calcBeta(resid,residAnt,z,SL->n);
      // printf("BETA %d: \n", beta);

      // // calcula prox direcao de busca
      // copiaMat(direc,dAnt,SL->n, 1); // pAnt = p
      // calcProxDirecBusca(direc,resid,beta,dAnt,SL->n); 
      // printf("p: \n");
      // prnMat(direc,SL->n,1);
    }

    return it;  
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

void prnMat (double **mat, unsigned int n, unsigned int m){
  for (unsigned int i=0; i < n;++i){
    for (unsigned int j=0; j < m;++j)
      printf(" %f",mat[i][j]);
    printf("\n");
  }
}

void copiaMat(double **matA, double **matB,int lin,int col){
  for (int i=0;i < lin;++i)
    for (int j=0;j < col;++j)
      matB[i][j] = matA[i][j]; 
}