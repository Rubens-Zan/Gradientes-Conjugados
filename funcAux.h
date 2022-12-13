#include <stdbool.h>

typedef struct tComando {
    int dimensao,nDiagonais, nIter, erroMax;
    bool usarPreCondicionador;   
    char saida[100];
} tComando;

void tratamentoEntrada(int argc, char **argv, tComando *comando);
