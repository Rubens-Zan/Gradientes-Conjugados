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
double generateRandomA(unsigned int i, unsigned int j, unsigned int k);
double generateRandomB(unsigned int k);