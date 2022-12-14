#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>
#include "utils.h"

// double timestamp(void)
// {
//   struct timespec tp;
//   clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
//   return((double)(tp.tv_sec + tp.tv_nsec*1.0e-9));
// }



/**
 * @brief Função para fazer o tratamento da entrada
 * @param argc 
 * @param argv 
 * @param comando 
 */
void tratamentoEntrada(int argc, char **argv, tComando *comando){
    for (int i=0;i<argc;++i){
        if (strcmp ( argv[i], "-n") == 0 && (i+1 < argc)){
            comando->dimensao = atoi(argv[i+1]); 
            i++; 
        }else if (strcmp ( argv[i], "-k") == 0 && (i+1 < argc)){
           comando->nDiagonais = atoi(argv[i+1]);
           i++;
        }else if (strcmp ( argv[i], "-p") == 0 && (i+1 < argc)){
            comando->usarPreCondicionador= argv[i+1] > 0; 
            i++;
        }else if(strcmp ( argv[i], "-i") == 0 && (i+1 < argc)){
            comando->nIter = atoi(argv[i+1]); 
            i++;
        }else if(strcmp ( argv[i], "-e") == 0 && (i+1 < argc)){
            comando->erroMax = atoi(argv[i+1]); 
            i++;           
        }else if (strcmp ( argv[i], "-o") == 0 && (i+1 < argc)){
           strcpy(comando->saida, argv[i+1]);
           i++;
        }  
    }

}