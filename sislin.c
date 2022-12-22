#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>

#include "opmatrizes.h"
#include "sislin.h"
#include "utils.h"
#include <string.h>

// Valor absoluto de um número. Alternativa ao uso da função 'fabs()'
#define ABS(num)  ((num) < 0.0 ? -(num) : (num))

// Alocaçao de matriz em memória.
SistLinear_t *alocaSisLin(unsigned int n)
{
  SistLinear_t *SL = (SistLinear_t *)malloc(sizeof(SistLinear_t));

  if (SL)
  {

    SL->n = n;
    SL->A = (real_t **)malloc(n * sizeof(real_t *));
    SL->b = (real_t *)malloc(n * sizeof(real_t));

    if (!(SL->A) || !(SL->b))
    {
      liberaSisLin(SL);
      return NULL;
    }

    // Matriz como vetor de N ponteiros para um único vetor com N*N elementos
    SL->A[0] = (real_t *)malloc(n * n * sizeof(real_t));
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
  unsigned int n = SL->n;
  // inicializando sequencia de numeros aleatorios
  srand(20222);

  // inicializa vetor b
  for (unsigned int i = 0; i < n; ++i)
  {
    SL->b[i] = generateRandomB(nDiagonais);
  }

  // inicializa a matriz A
  for (unsigned int i = 0; i < n; ++i)
  {
    for (unsigned int j = 0; j < n; ++j)
    {
      SL->A[i][j] = generateRandomA(i, j, nDiagonais);
    }
  }

  // FAZER SER DIAGONAL DOMINANTE
  //  FAZER O TERMO I,I SER A SOMA DE TUDO
  for (unsigned int i = 0; i < n; ++i)
  {
    SL->A[i][i] = 0.1;
    for (unsigned int j = 0; j < n; ++j)
    {
      SL->A[i][i] += SL->A[i][j];
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
void calcProxDirecBusca(double *proxDir, double *z, double beta, double *direcAnterior, int n)
{
  for (int i = 0; i < n; ++i)
  {
    proxDir[i] = z[i] + (beta * direcAnterior[i]);
  }
}

/**
 * @brief - Calcula alpha
 *
 * alpha = r<k>^t * z<k> / (p<k>^t * A * p<k>)
 * ok
 * @param resid - residuo calculado
 * @param A - Matriz original
 * @param p - Matriz de direção de busca calculada
 * @param z - Matriz z
 * @param n - dimensão da matriz
 * @return double
 */
double calcAlpha(double *resid, double **A, double *p, double *z, int n)
{
  double alpha = 0;
  double *vetorAux = (double *)malloc(sizeof(double) * n);

  double resTxZ = 0;

  for (int i = 0; i < n; ++i)
  { // resTxZ = resid^T * z
    resTxZ += resid[i] * z[i];
  }

  for (int i = 0; i < n; ++i)
  { // multPtxA = pT * A
    vetorAux[i] = 0;
    for (int j = 0; j < n; ++j)
    {
      vetorAux[i] = p[j] * A[j][i];
    }
  }

  // multiplicaMatriz_Vetor(A,p,vetorAux,n ); // multPtxA = p * A

  double denom = 0;
  //  multiplicaVetor_Vetor(vetorAux,p,n);
  for (int i = 0; i < n; ++i)
  {
    denom += vetorAux[i] * p[i];
  }

  alpha = resTxZ / denom;
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
 * @param resid
 * @param residAnt
 * @param z
 * @param zAnt
 * @param n
 * @return double
 */
double calcBeta(double *resid, double *residAnt, double *z, double *zAnt, int n)
{
  double beta = 0;
  double residTxZ = 0;
  for (int i = 0; i < n; ++i)
  { // resTxZ = resid^T * z
    residTxZ += resid[i] * z[i];
  }
  double resTantxZ = 0;
  for (int i = 0; i < n; ++i)
  { // res^Tant * zAnt
    resTantxZ += residAnt[i] * zAnt[i];
  }

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
void calcX(double *proxX, double *xAnt, double alpha, double *p, int n)
{
  for (int i = 0; i < n; ++i)
  {
    proxX[i] = xAnt[i] + (alpha * p[i]);
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
void calcResiduo(double *residuoAnterior, double alpha, double **A, double *p, double *residuo, int n)
{
  double *vetorAux = (double *)malloc(sizeof(double) * n);

  // alpha * A * p<k>
  for (int i = 0; i < n; ++i)
  { // OK
    vetorAux[i] = 0;
    for (int j = 0; j < n; j++)
    {
      vetorAux[i] += (alpha * A[i][j]) * p[j];
    }
  }

  for (int i = 0; i < n; ++i)
  {
    residuo[i] = residuoAnterior[i] - vetorAux[i]; // ok
  }
}

void calcZ(real_t *z, real_t *inverse_c, real_t *residuo, unsigned int size)
{
  for (int i = 0; i < size; ++i)
  {
    z[i] = inverse_c[i] * residuo[i];
    // veriifcar nan
  }
}

/**
 * @brief Função para formar Sistema linear para ficar no formato (A^T) * A * x = (A^T) * b, facilita a conversão
 *
 * @param SL
 */
void formataSLGradConj(SistLinear_t *SL)
{
  // A = A^T * A
  // b = A^T * b
  double **matAux = alocarMatriz(SL->n, SL->n);
  transporMat(SL->A, matAux, SL->n, SL->n);
  multMat(SL->A, matAux, SL->n, SL->n, SL->n, SL->n, SL->A);
  multiplicaMatriz_Vetor(matAux, SL->b, SL->b, SL->n); // poder estar errado
};

//////////////////////////////////////
void precondicionador_jacobi(SistLinear_t *SL, real_t *M)
{
  for (int i = 0; i < SL->n; ++i)
  {
    M[i] = (real_t)1 / SL->A[i][i]; // gera vetor com diagonais do SL
  }
}

void precondicionador_identidade(SistLinear_t *SL, real_t *M)
{
  for (int i = 0; i < SL->n; i++)
  {
    M[i] = 1;
  }
}

void aplicaPreCondicSL(SistLinear_t *SL, real_t *M)
{
  for (int i = 0; i < SL->n; i++)
  {
    for (int j = 0; j < SL->n; j++)
    {
      if (i == j)
      {
        SL->A[i][j] = SL->A[i][j] * M[i]; // aplica pre condicionador se for na diagonal
      }
      if (i == j)
      {
        SL->b[i] *= M[i];
      }
    }
  }
}

/******************NORMAS**********************************/
/**
 * @brief 
 *  r = b - Ax
 * @param b 
 * @param A 
 * @param x 
 * @return double 
 */
double calcErroNormaEuc(double *b, double **A, double *x, int n){
  double *vAux = (double *)malloc(sizeof(double) * n);     // matriz de chute anterior
  double normaEucl = 0;
  multiplicaMatriz_Vetor(A,x,vAux, n);
  subtraiVetor(b,vAux, vAux,n);
  normaEucl = sqrt(multiplicaVetor_Vetor(vAux, vAux, n));

  return normaEucl;
}

double calcNormaMaxRel(double *xAnt,double *x, int n){
  double normaMaxRel = 0; 
  for (int i=0;i < n;++i){
    if ( (x[i] - xAnt[i]/ x[i]) > normaMaxRel){
      normaMaxRel =  ABS(x[i] - xAnt[i])/ABS(x[i]); 
    }
  }

  return normaMaxRel; 
}

/******************FUNCOES GRADIENTE CONJUGADO E PRE CONJUGADO**********************************/

/**
 * @brief Metodo resolvedor de sistemas lineares pelo metodo de gradiente conjugado com uso do pre condicionador de jacobi
 *
 * @param SL
 * @param maxIt
 * @param tol
 * @param matSaida
 * @param arqSaida
 * @return int
 */
int gradienteConjugadoPreCondic(SistLinear_t *SL, int maxIt, double tol, double matSaida[][2], FILE *arqSaida)
{
  // inicia chute inicial com vetor de 0
  double alpha, beta;
  double *resid = (double *)malloc(sizeof(double) * SL->n);    // matriz de residuo
  double *residAnt = (double *)malloc(sizeof(double) * SL->n); // matriz de residuo anterior
  double *direc = (double *)malloc(sizeof(double) * SL->n);    // matriz de direcao de busca
  double *direcAnt = (double *)malloc(sizeof(double) * SL->n); // matriz de direcao anterior
  double *xAnt = (double *)malloc(sizeof(double) * SL->n);     // matriz de chute anterior
  double *z = (double *)malloc(sizeof(double) * SL->n);        // matriz de
  double *zAnt = (double *)malloc(sizeof(double) * SL->n);     // matriz de
  double *x = (double *)malloc(sizeof(double) * SL->n);

  double *auxMatJacobi = (double *)malloc(sizeof(double) * SL->n);

  double tMedioIter, tempoResid, tempoPreCond;
  tMedioIter = 0;
  tempoResid = timestamp();
  ;

  int it;

  // Ajusta tipos iniciais
  // A=(A^T) * A
  // b=(A^T) * b
  formataSLGradConj(SL);
  // (M^-1) * A
  // (M^-1) * b
  precondicionador_jacobi(SL, auxMatJacobi);
  aplicaPreCondicSL(SL, auxMatJacobi);
  tempoPreCond = timestamp() - tempoResid;
  // residuo = b
  copiaVetor(SL->b, resid, SL->n); // como x inicial igual a 0, desconsidero o r = b - (A * x)
  // calcula z
  calcZ(z, auxMatJacobi, resid, SL->n);
  // direc = z
  copiaVetor(z, direc, SL->n);

  // verifica erro do x com a tolerancia
  // loop
  for (it = 0; it < maxIt; ++it)
  {
    double tIterInicio = timestamp();
    

    // ALPHA = z * r / pT * A * p
    alpha = calcAlpha(resid, SL->A, direc, z, SL->n);
    // calcula novo x
    copiaVetor(x, xAnt, SL->n); // xant = x
    calcX(x, xAnt, alpha, direc, SL->n);

    double normaMaxRel =calcNormaMaxRel(xAnt, x, SL->n); 
    fprintf(arqSaida, "# iter %d: %.15g", it, normaMaxRel);
    // calcula novo residuo
    copiaVetor(resid, residAnt, SL->n);
    calcResiduo(residAnt, alpha, SL->A, direc, resid, SL->n);

    // if erro
    // o erro aproximado em x após a k-ésima iteração (max|xi - xi-1|);
    if (tol != 0 && tol > normaMaxRel)
    {
      tMedioIter += timestamp() - tIterInicio;
      break;
    }
    // z0 = z
    copiaVetor(z, zAnt, SL->n);
    // calcula z
    calcZ(z, auxMatJacobi, resid, SL->n);
    // beta = rt * z / rtant *zant
    beta = calcBeta(resid, residAnt, z, zAnt, SL->n);
    // calcula prox direcao
    copiaVetor(direc, direcAnt, SL->n); // dAnt = d
    calcProxDirecBusca(direc, z, beta, direcAnt, SL->n);

    tMedioIter += timestamp() - tIterInicio;
  }

  // A norma euclidiana do resíduo (||r||), onde r = b - Ax
  fprintf(arqSaida, "# residuo: || %.15g ||", calcErroNormaEuc(SL->b,SL->A, x, SL->n));

  fprintf(arqSaida, "# Tempo PC: %.15g \n", tempoPreCond);
  tMedioIter = tMedioIter / it;
  fprintf(arqSaida, "# Tempo tMedioIter:: %.15g \n", tMedioIter);
  tempoResid = timestamp() - tempoResid;
  fprintf(arqSaida, "# Tempo residuo:: %.15g", tempoResid);

  free(resid);
  free(residAnt);
  free(direc);
  free(direcAnt);
  free(xAnt);
  free(x);
  free(z);
  free(zAnt);

  return it;
}

/**
 * @brief Metodo resolvedor de sistemas lineares pelo metodo de gradiente conjugado
 *
 * @param SL
 * @param maxIt
 * @param tol
 * @param matSaida
 * @param arqSaida
 * @return int
 */
int gradienteConjugado(SistLinear_t *SL, int maxIt, double tol, double matSaida[][2], FILE *arqSaida)
{
  double tMedioIter, tempoResid, tempoPreCond;
  double alpha, beta;
  double *resid = (double *)malloc(sizeof(double) * SL->n);    // matriz de residuo
  double *residAnt = (double *)malloc(sizeof(double) * SL->n); // matriz de residuo anterior
  double *direc = (double *)malloc(sizeof(double) * SL->n);    // matriz de direcao de busca
  double *direcAnt = (double *)malloc(sizeof(double) * SL->n); // matriz de direcao anterior
  double *xAnt = (double *)malloc(sizeof(double) * SL->n);     // matriz de chute anterior
  double *x = (double *)malloc(sizeof(double) * SL->n);
  double *vetorAux = (double *)malloc(sizeof(double) * SL->n);
  int it = 0;
  tMedioIter = 0;

  tempoResid = timestamp();
  // transforma SL em ax=b em (A^T)Ax = A^t *b
  formataSLGradConj(SL);
  
  tempoPreCond = timestamp() - tempoResid;

  // residuo = b
  copiaVetor(SL->b, resid, SL->n);
  // direcao = residuo
  copiaVetor(resid, direc, SL->n);

  // for
  for (it = 0; it < maxIt; ++it)
  {
    double tIterInicio = timestamp();

    // CALCULA APLHA
    // Alpha = (residuo^t * residuo) / (p^t * A * p )
    for (int i = 0; i < SL->n; ++i)
    { // multPtxA = pT * A
      vetorAux[i] = 0;
      for (int j = 0; j < SL->n; ++j)
      {
        vetorAux[i] = direc[j] * SL->A[j][i];
      }
    }
    alpha = (multiplicaVetor_Vetor(resid, resid, SL->n)) / (multiplicaVetor_Vetor(vetorAux, direc, SL->n));
    // calcula novo x
    // x1 = x0 +alpha * p
    copiaVetor(x, xAnt, SL->n);
    calcX(x, xAnt, alpha, direc, SL->n);
    double normaMaxRel =calcNormaMaxRel(xAnt, x, SL->n); 
    fprintf(arqSaida, "# iter %d: %.15g", it, normaMaxRel);

    // calcular novo residuo
    // r1= r0 - alpha0 * A *p0
    copiaVetor(resid, residAnt, SL->n);
    calcResiduo(residAnt, alpha, SL->A, direc, resid, SL->n);

    // normaMaxErroRelativo
    // if
    if (tol != 0 && tol > normaMaxRel)
    {
      tMedioIter += timestamp() - tIterInicio;
      break;
    }
  
    // escalarBeta = vetor residuo^t * vetor residuo / vetor residuo anteriorT * vetor residuo anterior
    beta = multiplicaVetor_Vetor(resid, resid, SL->n) / multiplicaVetor_Vetor(residAnt, residAnt, SL->n);
  
    // calcula prox direc
    // P = R + b0 * pAnt
    copiaVetor(direc, direcAnt, SL->n);
    calcProxDirecBusca(direc, resid, beta, direcAnt, SL->n);

    tMedioIter += timestamp() - tIterInicio;
  }

  fprintf(arqSaida, "# residuo: || %.15g ||", calcErroNormaEuc(SL->b,SL->A, x, SL->n));
  
  // tempo final
  fprintf(arqSaida, "# Tempo PC: %.15g \n", tempoPreCond);
  tMedioIter = tMedioIter / it;
  fprintf(arqSaida, "# Tempo tMedioIter:: %.15g \n", tMedioIter);
  tempoResid = timestamp() - tempoResid;
  fprintf(arqSaida, "# Tempo residuo:: %.15g", tempoResid);
  free(resid);
  free(residAnt);
  free(direc);
  free(direcAnt);
  free(xAnt);
  free(x);

  return it;
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

void prnVetor(real_t *v, unsigned int n)
{
  int i;

  printf("\n");
  for (i = 0; i < n; ++i)
    printf("%10g ", v[i]);
  printf("\n\n");
}

void prnVetorArq(real_t *v, unsigned int n, FILE *arqSaida)
{
  int i;
  fprintf(arqSaida, "\n");
  for (i = 0; i < n; ++i)
  {
    fprintf(arqSaida, "%.15g \n", v[i]);
  }
  fprintf(arqSaida, "\n\n");
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
