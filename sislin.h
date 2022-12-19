#ifndef __SISLIN_H__
#define __SISLIN_H__

#define COEF_MAX 32.0 // Valor máximo usado para gerar valores aleatórios de
		      // coeficientes nos sistemas lineares.

// real_t: tipo usado para representar valores em ponto flutuante
typedef double real_t;

// rtime_t: tipo usado para representar valores de tempo em ponto flutuante
typedef double rtime_t;

// Estrutura para definiçao de um sistema linear qualquer
typedef struct {
  real_t **A; // coeficientes
  real_t *b; // termos independentes
  unsigned int n; // tamanho do SL
} SistLinear_t;


// Alocaçao e desalocação de matrizes
SistLinear_t* alocaSisLin (unsigned int n);
void liberaSisLin (SistLinear_t *SL);
void iniSisLin (SistLinear_t *SL, unsigned int nDiagonais);

//
double calcBeta(double **resid,double **residAnt, int n); 
double calcAlpha(double **resid,double **A, double **p, int n);
void calcProxX(double **proxX,double **xAnt,double alpha, double **p, int n); 
void calcProxDirecBusca(double **proxDir,double **resid, double beta,double **direcAnterior, int n); 
void calcResiduo(double **residuoAnterior, double alpha, double **A, double **p, double ** residuo,int n);
//

void prnMat (double **mat, unsigned int n, unsigned int m);
// Leitura e impressão de sistemas lineares
SistLinear_t *lerSisLin ();
void prnSisLin (SistLinear_t *SL);
void prnVetor (real_t *vet, unsigned int n);
int gradienteConjugado(SistLinear_t *SL, double **x, double *M, int maxIt, double tol, double matSaida[][2]);
void copiaMat(double **matA, double **matB,int lin,int col); 
#endif // __SISLIN_H__

