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
 * p[i][j]<k> = r[i][j]<k> + beta<k-1> * p[i][j]<k-1> 
 * @param proxDir
 * @param resid - Residuo calculado r<i>
 * @param beta - Beta calculado beta<i-1>
 * @param direcAnterior - Direcao anterior calculada 
 * @return double 
 */
void calcProxDirecBusca(double **proxDir,double **resid, double beta,double **direcAnterior, int n){
    for (int i=0; i < n;++i)
      for(int j=0;j < 1;++j)
        proxDir[i][j] = resid[i][j] + beta * direcAnterior[i][j]; 
}

/**
 * @brief 
 * alpha = r<k>^t*r<k>/(p<k>^t*A*p<k>)
 * @param resid - residuo calculado
 * @param A - Matriz original
 * @param p - Matriz de direção de busca calculada
 * @param n - dimensão da matriz
 * @return double 
 */
double calcAlpha(double **resid,double **A, double **p, int n){
  double alpha = 0; 
  double **residTransp=alocarMatriz(1,n+1); // matriz de residuo transposta
  double **pTransp = alocarMatriz(1,n+1); // matriz de direcao transposta
  transporMat(resid, residTransp, n, 1);
  transporMat(p, pTransp, n, 1);
  
  double **resMultResult = multMat(resid,residTransp,n,1,1,n); // resid * resid^T 
  double **multResulPtxA = multMat(pTransp,A,1,n,n,n); // p^T * A
  double  **multResulPxAxP =  multMat(multResulPtxA,p,1,n,n,1); // (p^T * A) * p

  alpha = resMultResult[0][0] /  multResulPxAxP[0][0]; 

  liberarMatriz(residTransp);
  liberarMatriz(pTransp);
  liberarMatriz(resMultResult);
  liberarMatriz(multResulPtxA);
  liberarMatriz(multResulPxAxP);
  
  return alpha; 
}

/**
 * @brief Calcula beta com base no residuo atual e no anterior
 * beta<k> = (residuo<k+1>^T*residuo<k+1>)/(residuoAnterior<k>^T*residuoAnterior<k>)
 * @param resid - Residuo da k-esima iteração
 * @param residAnt - Resíduo da k-1-esima iteração
 * @return double - valor de beta calculado
 */
double calcBeta(double **resid,double **residAnt, int n){
	double beta = 0; 
  double **residTransp=alocarMatriz(1,n+1); // matriz de residuo transposta
	double **residAntTransp=alocarMatriz(1,n+1); // matriz do residuo anterior transposta
  
  transporMat(resid, residTransp, n,1);
  transporMat(residAnt, residAntTransp, n,1);
  double **resMultResult = multMat(resid,residTransp,n,1,1,n); // res * res^t
  double **resAntMultResult = multMat(residAnt,residAntTransp,n,1,1,n); // res<k-1> * res<k-1>^T

  beta += resMultResult[0][0]/ resAntMultResult[0][0];

  liberarMatriz(resMultResult);
  liberarMatriz(resAntMultResult);
  liberarMatriz(residAntTransp);
  liberarMatriz(residTransp);

  return beta; 
}


/**
 * @brief Calcula o proximo x<k>
 * x[i][j]<k> = x[i][j]<k-1> + a <k> * p[i][j]<k-1> 
 * @param proxX  
 * @param xAnt - Valor do x<i-1> anterior
 * @param alpha - valor do alpha<i-1> anterior
 * @param p - valor da direção de busca p<i-1)> anterior
 * @param n - dimensao do sistema linear
 * @return double - Proximo x calculado
 */
void calcProxX(double **proxX,double **xAnt,double alpha, double **p, int n){
    for (int i =0;i < n;++i)
      for (int j=0;j < 1;++j)
        proxX[i][j] = xAnt[i][j] + alpha * p[i][j]; 
}

/**
 * @brief 
 * r<k> = r<k-1> -alpha<k-1> * A * p<k-1>
 * @param residuoAnterior 
 * @param alpha 
 * @param A 
 * @param p 
 * @param residuo 
 * @param n - dimensao do sistema linear
 */
void calcResiduo(double **residuoAnterior, double alpha, double **A, double **p, double ** residuo,int n){
  double **matAlphaxA = alocarMatriz(n+1,n+1);

  for (int i=0;i < n;++i){
    for (int j=0;j<n;++j){
      matAlphaxA[i][j] = alpha * A[i][j]; // alpha * A
    }
  }
  double **multAlphaxAxp =multMat(matAlphaxA, p,n,n,n,1); // (alpha * A) * p  

  for (int i=0;i < n;++i)
    for(int j=0;j < 1;++j){
      residuo[i][j] = residuoAnterior[i][j]-multAlphaxAxp[i][j];  
    }

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
  double **matAux = multMat(A, x, n,n,n,1);
  
  for (int i=0;i < n;++i){
    for (int j=0;j < 1; ++j){
      resid[i][j] =b[i] - matAux[i][j];
    }
  }
}

/**
 * @brief Metodo resolvedor de sistemas lineares pelo metodo de gradiente conjugado
 * 
 * @param A  - MATRIZ A SER TESTADA
 * @param b - Coeficiente de solucao
 * @param x -  Solucao do sistema inicial
 * @param M  - Matriz precondicionadora 
 * @param maxIt - NUMERO MAXIMO DE ITERACOES
 * @param tol - TOLERANCIA
 * @param n  - DIMENSAO 
 * @param matSaida - Matriz que guardara o Erro e a Norma de cada iteração
 * @return int - NUMERO DE ITERACOES
 */
int gradienteConjugado(SistLinear_t *SL, double **x, double *M, int maxIt, double tol, double matSaida[][2]){
    // inicia chute inicial com vetor de 0  
    double alpha, beta; 
    double **resid = alocarMatriz(SL->n+1,2); // matriz de residuo 
    double **residAnt = alocarMatriz(SL->n+1,2); // matriz de direcao
    double **p = alocarMatriz(SL->n+1,2); // matriz de direcao
    double **pAnt = alocarMatriz(SL->n+1,2); // matriz de direcao    
    double **xAnt = alocarMatriz(SL->n+1,2); // matriz de direcao
    double **z = alocarMatriz(SL->n+1,2); // matriz de direcao
  
    int it;
    
    calcResiduoInicial(SL->A,SL->b,x,resid,SL->n); 
    
    // inicia direcao inicial com o residuo inicial
    copiaMat(resid,p,SL->n, 1); 

    // verifica erro do x com a tolerancia 
    // loop 
    for(it =0;it < maxIt;++it){
      // calcula alpha
      alpha = calcAlpha(resid,SL->A,p,SL->n);
      // printf("ALPHA %f\n", alpha);
      copiaMat(x,xAnt,SL->n,1); // xAnterior = xAtual
      // calcula x
      calcProxX(x,xAnt,alpha,p,SL->n); 
      // printf("X: \n");
      // prnMat(x,SL->n,1);

      // calcula residuo
      copiaMat(resid,residAnt,SL->n, 1); 
      calcResiduo(residAnt,alpha,SL->A,p,resid,SL->n); 
      printf("residuo: \n");
      prnMat(resid,SL->n,1);
      
      // calcula beta
      beta = calcBeta(resid,residAnt,SL->n);
      // printf("beta %f\n", beta);

      // calcula prox direcao de busca
      copiaMat(p,pAnt,SL->n, 1); // pAnt = p
      calcProxDirecBusca(p,resid,beta,pAnt,SL->n); 
      // printf("p: \n");
      // prnMat(p,SL->n,1);
    }

    return it;  
}


/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
  \param residue Valor do resíduo

  \return Norma L2 do resíduo
*/
// real_t normaL2Residuo(SistLinear_t *SL, real_t *x)
// {
//   real_t norma = 0;
//   real_t res[SL->n + 1];
//   // A*X
//   real_t matrix_solution[SL->n + 1];

//   // função que calcula a Multiplicação A*X da matriz
//  real_t ** matrix_solution = multMat(SL, x);

//   // r0 = b0 - x0, r1 = b1 - x1, r2 = b2 - x2.... rn = bn - xn
//   for (int i = 0; i < SL->n; i++)
//     res[i] = SL->b[i] - matrix_solution[i];

//   // norma = sqrt(r1²+r2²+r3²+...)
//   for (int i = 0; i < SL->n; i++)
//     norma += res[i] * res[i];
//   norma = sqrt(norma);
//   // Retorna número de iterações
//   if (isnan(norma) || isinf(norma))
//   {
//     fprintf(stderr, "Resultado invalido! \n");

//     return -1;
//   }
//   return norma;
// }



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