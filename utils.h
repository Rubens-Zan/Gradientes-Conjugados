#include <stdbool.h>

typedef struct tComando {
    int dimensao,nDiagonais, nIter;
    double erroMax;
    bool usarPreCondicionador;   
    char saida[100];
} tComando;
double calcNormaMaxRel(double *xAnt,double *x, int n);
double normaMaxErroRelativo(double *x, double *xAnt, unsigned int n);

double calcularNormaL2Residuo( double *residuo, unsigned int n);

double timestamp(void);
void tratamentoEntrada(int argc, char **argv, tComando *comando);
double ** alocarMatriz(int lin,int col); 
void liberarMatriz(double **matriz);
void inicializarMatriz(double **vet, int lin,int col);

