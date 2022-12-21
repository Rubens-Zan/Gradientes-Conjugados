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
 * p<k> = z<k> + beta<k-1> * p<k-1>  
 * 
 * @param proxDir
 * @param z - z calculado 
 * @param beta - Beta calculado beta<i-1>
 * @param direcAnterior - Direcao anterior calculada 
 * @return double 
 */
void calcProxDirecBusca(double *proxDir,double *z, double beta,double *direcAnterior, int n){
  for (int i=0; i < n;++i){
    proxDir[i] = z[i] + (beta * direcAnterior[i]); 
  }
}

/**
 * @brief - Calcula alpha
 * 
 * alpha = r<k>^t * z<k> / (p<k>^t * A * p<k>)  
 * 
 * @param resid - residuo calculado
 * @param A - Matriz original
 * @param p - Matriz de direção de busca calculada
 * @param z - Matriz z
 * @param n - dimensão da matriz
 * @return double 
 */
double calcAlpha(double *resid,double **A, double *p,double *z, int n){
  double alpha = 0; 
  double *vetorAux = (double *) malloc (sizeof(double)* n);

  double resTxZ = multiplicaVetor_Vetor(resid,z,n); // resTxZ = resid^T * z

  multiplicaMatriz_Vetor(A,p,vetorAux,n ); // multPtxA = p * A
  double denom =multiplicaVetor_Vetor(vetorAux,p,n); 

  alpha = resTxZ /  denom; 
  // Verificação se resultou em NaN ou +/- infinito
  if (isnan(alpha) || isinf(alpha))
  {
    fprintf(stderr, "Erro alpha(calcAlpha): %g é NaN ou +/-Infinito\n", alpha);
    exit(1);
  }
  free(vetorAux);

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
double calcBeta(double *resid,double *residAnt,double *z, int n){
	double beta = 0; 
  double residTxZ = multiplicaVetor_Vetor(resid, z, n); // residuo * z
  double resTantxZ = multiplicaVetor_Vetor(residAnt,z,n); // residuoAnterior * z 

  beta = residTxZ / resTantxZ;
  // Verificação se resultou em NaN ou +/- infinito
  if (isnan(beta) || isinf(beta))
    {
      fprintf(stderr, "Erro beta(calcbeta): %g é NaN ou +/-Infinito\n", beta);
      exit(1);
    }
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
void calcX(double *proxX,double *xAnt,double alpha, double *p, int n){
  double *vetorAux = (double *) malloc (sizeof(double)* n);
  multiplicaInteiro_Vetor(alpha, p,vetorAux,n); //vet aux =  alpha * p
  somaVetor(xAnt,vetorAux,proxX, n); // proxX = vet aux + x ant
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
void calcResiduo(double *residuoAnterior, double alpha, double **A, double *p, double *residuo,int n){
  double *vetorAux = (double *) malloc (sizeof(double)* n);
  multiplicaMatriz_Vetor(A,p,vetorAux,n);
  multiplicaInteiro_Vetor(alpha,vetorAux, vetorAux, n);
  subtraiVetor(residuoAnterior,vetorAux, residuo,n); 
}

void calcZ(real_t *z, real_t *inverse_c, real_t *residuo,unsigned int size){
  for(int i=0;i<size;++i){
    z[i] = inverse_c[i] * residuo[i];
    // veriifcar nan 
  }
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
int gradienteConjugadoPreCondic(SistLinear_t *SL, double *x, double *matPreConj, int maxIt, double tol, double matSaida[][2]){
    // inicia chute inicial com vetor de 0  
    double alpha, beta; 
    double *resid = (double *) malloc (sizeof(double)*SL->n); // matriz de residuo 
    double *residAnt = (double *) malloc (sizeof(double)*SL->n); // matriz de residuo anterior
    double *direc = (double *) malloc (sizeof(double)*SL->n); // matriz de direcao de busca
    double *dAnt = (double *) malloc (sizeof(double)*SL->n); // matriz de direcao anterior    
    double *xAnt = (double *) malloc (sizeof(double)*SL->n); // matriz de chute anterior
    double *z = (double *) malloc (sizeof(double)*SL->n); // matriz de 
  
    int it;
    // as inicializacoes estao ok!    
    // x = 0 
    // inicializarMatriz(x,SL->n,1);  
    memset(x,0,sizeof(x)*SL->n);

    // r<0> = b - A * x<0>
    // calcResiduoInicial(SL->A,SL->b,x,resid,SL->n);
    copiaVetor(SL->b,resid,SL->n); // como chute inicial =0, copia o vetor b 

    // z<0> = MATPRECONJ^-1 * r<0>
    calcZ(z,matPreConj,resid,SL->n); 

    // inicia direcao inicial com o z inicial]
    copiaVetor(z, direc,SL->n); 

    // verifica erro do x com a tolerancia 
    // loop 
    for(it =0;it < maxIt;++it){
      // calcula alpha
      alpha = calcAlpha(resid,SL->A,direc,z,SL->n);

      copiaVetor(x, xAnt,SL->n); 
      
      // calcula x
      calcX(x,xAnt,alpha, direc, SL->n);

      // calcula residuo
      copiaVetor(resid, residAnt,SL->n); 

      calcResiduo(residAnt,alpha,SL->A,direc,resid,SL->n); 
      
      // // calcula z
      // z<k+1> = C^-1 * r<k+1>
      calcZ(z, matPreConj,resid, SL->n); 
      
      // // calcula erro 
      // // ERRO = r<k+1> * r<k+1>

      // // faz o if com o  erro e tolerancia 
      
      // // calcula beta
      beta = calcBeta(resid,residAnt,z,SL->n);

      // // calcula prox direcao de busca
      copiaVetor(direc, dAnt,SL->n); 
      
      calcProxDirecBusca(direc,z,beta,dAnt,SL->n); 
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


//////////////////////////////////////
void precondicionador_jacobi(SistLinear_t *SL, real_t *M){
  for (int i=0;i < SL->n;++i){
    M[i]= (real_t) 1 / SL->A[i][i]; // gera vetor com diagonais do SL
  }
}

void precondicionador_identidade(SistLinear_t *SL, real_t *M){
  for(int i = 0; i < SL->n;i++) {
    M[i] = 1;
	}
}


void aplicaPreCondicMat(SistLinear_t *SL, real_t *M){
  for (int i = 0; i < SL->n; i++){
    for (int j = 0; j < SL->n; j++){
      if (i == j)
        SL->A[i][j] = SL->A[i][j] * M[i];  //aplica pre condicionador se for na diagonal
    }
  }
}