#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>

#define COEF_MAX 32.0 // Valor máximo usado para gerar valores aleatórios de
		      // coeficientes nos sistemas lineares.


// Estrutura para definiçao de um sistema linear qualquer
typedef struct {
  double **A; // coeficientes
  double *b; // termos independentes
  unsigned int n; // tamanho do SL
} SistLinear_t;


// Alocaçao e desalocação de matrizes
SistLinear_t* alocaSisLin (unsigned int n);
void liberaSisLin (SistLinear_t *SL);
void iniSisLin (SistLinear_t *SL, unsigned int nDiagonais);

//

//
void prnMat (double **mat, unsigned int n, unsigned int m);
// Leitura e impressão de sistemas lineares
void prnSisLin (SistLinear_t *SL);
void prnVetor (double *vet, unsigned int n);
double multiplicaVetores(double *vetA, double *vetB, unsigned int n);
void copiaVetor (double *a, double *b, unsigned int N);
#endif // __SISLIN_H__

