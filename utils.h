#include <stdbool.h>

typedef struct tComando {
    int dimensao,nDiagonais, nIter, erroMax;
    bool usarPreCondicionador;   
    char saida[100];
} tComando;

void tratamentoEntrada(int argc, char **argv, tComando *comando);
double ** alocarMatriz(int lin,int col); 
void liberarMatriz(double **matriz);
void inicializarMatriz(double **vet, int lin,int col);

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
