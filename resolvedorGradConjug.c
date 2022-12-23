#include "sislin.h"
#include "utils.h"
#include "resolvedorGradConjug.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief - Calcula alpha parece ok 
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
double calcAlphaPreCond(double *resid, double **A, double *p, double *z, int n)
{
    double alpha = 0;
    double *pTxA = (double *)malloc(sizeof(double) * n); // vetor aux para calculo de p^T * A

    double residTxZ = multiplicaVetores(z, resid, n);

    //Percorre a matriz SL->A e o vetor D
    for (int i = 0; i < n; ++i)
    {
        double soma = 0.0;
        for (int j = 0; j < n; ++j)
        {
            //calcula o iésimo elemento do vetor_resultante (esse vetor é chamado de 'z' no livro M.Cristina C. Cunha)
            soma = soma + A[i][j] * z[j];

            // Teste para ver se não foi gerado um NaN ou um número infinito
            if (isnan(soma) || isinf(soma))
            {
                fprintf(stderr, "Erro soma(calcAlphaPreCond): %g é NaN ou +/-Infinito\n", soma);
                exit(1);
            }
        }
        // O iésimo elemento do vetor resultante recebe soma.
        pTxA[i] = soma;
    }

    // pT * A * p
    double pTxAxP = 0.0;
    for (int i = 0; i < n; ++i)
    {
        pTxAxP += z[i] * pTxA[i];
        // Teste para ver se não foi gerado um NaN ou um número infinito
        if (isnan(pTxAxP) || isinf(pTxAxP)){
            fprintf(stderr, "Erro pTxAxP(calcularDenominadorEscalarA): %g é NaN ou +/-Infinito\n", pTxAxP);
            exit(1);
        }
    }

    alpha = residTxZ / pTxAxP;
    // Verificação se resultou em NaN ou +/- infinito
    if (isnan(alpha) || isinf(alpha))
    {
        fprintf(stderr, "Erro alpha(calcAlphaPreCond): %g é NaN ou +/-Infinito\n", alpha);
        exit(1);
    }
    free(pTxA);

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
double calcBetaPreCond(double *resid, double *residAnt, double *z, double *zAnt, int n)
{
    double beta = 0;
    double residTxZ = multiplicaVetores(z,resid,n);
    double zAntxResidAnt = multiplicaVetores(zAnt, residAnt, n);  

    beta = residTxZ / zAntxResidAnt;
    // Verificação se resultou em NaN ou +/- infinito
    if (isnan(beta) || isinf(beta))
    {
        fprintf(stderr, "Erro beta(calcbetaPreCond): %g é NaN ou +/-Infinito\n", beta);
        exit(1);
    }
    return beta;
}

/**
 * @brief
 *
 * @param z
 * @param inverse_c
 * @param residuo
 * @param size
 */
void calcZ(double *z, double *inverse_c, double *residuo, unsigned int size)
{
    for (int i = 0; i < size; ++i)
    {
        z[i] = inverse_c[i] * residuo[i];
        // veriifcar nan
    }
}

//////////////////////////////////////

/**
 * @brief Inicializa a matriz M com a diagonal principal invertida de A
 *
 * @param SL
 * @param M
 */
void inicializaPreCondJacobi(SistLinear_t *SL, double *M)
{
    for (int i = 0; i < SL->n; ++i)
    {
        M[i] = (double)1 / SL->A[i][i]; // gera vetor com diagonais do SL
    }
}

void aplicaPreCondicSL(SistLinear_t *SL, double *M)
{
    for (int i = 0; i < SL->n; i++)
    {
        for (int j = 0; j < SL->n; j++)
        {
            SL->A[i][j] *= M[i]; // aplica pre condicionador se for na diagonal
            SL->b[i] *= M[i];
        }
    }
}

void inicializaSol(double *x, unsigned int n){
    for (int i=0;i < n;++i)
        x[i]=0;
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

    int it;

    inicializaPreCondJacobi(SL, auxMatJacobi);
 
    // (M^-1) * A
    // (M^-1) * b
    aplicaPreCondicSL(SL, auxMatJacobi);
    tempoPreCond = timestamp() - tempoResid;
    inicializaSol(x,SL->n);
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
        alpha = calcAlphaPreCond(resid, SL->A, direc, z, SL->n);
        // calcula novo x
        copiaVetor(x, xAnt, SL->n); // xant = x
        calcX(x, xAnt, alpha, direc, SL->n);

        double normaMaxRel = normaMaxErroRelativo(x,xAnt, SL->n);
        fprintf(arqSaida, "# iter %d: %.15g\n", it, normaMaxRel);
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
        copiaVetor(z, zAnt, SL->n); //////// verificar com z0 = r
        // calcula z
        calcZ(z, auxMatJacobi, resid, SL->n);
        // beta = rt * z / rtant *zant
        beta = calcBetaPreCond(resid, residAnt, z, zAnt, SL->n);
        // calcula prox direcao
        copiaVetor(direc, direcAnt, SL->n); // dAnt = d
        calcProxDirecBusca(direc, z, beta, direcAnt, SL->n);

        tMedioIter += timestamp() - tIterInicio;
    }

    // A norma euclidiana do resíduo (||r||), onde r = b - Ax
    fprintf(arqSaida, "# residuo: || %.15g || \n", calcularNormaL2Residuo(resid, SL->n));

    fprintf(arqSaida, "# Tempo PC: %.15g \n", tempoPreCond);
    tMedioIter = tMedioIter / it;
    fprintf(arqSaida, "# Tempo tMedioIter:: %.15g \n", tMedioIter);
    tempoResid = timestamp() - tempoResid;
    fprintf(arqSaida, "# Tempo residuo:: %.15g \n", tempoResid);

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

////////////////FUNCOES AUXILIAR SEM PRE CONDICIONADOR///////////////////////////

/**
 * @brief  OK
 * alpha = rT*r / direcT * A * direc
 * @param resid
 * @param direc
 * @param SL
 * @return double
 */
double calcAlpha(double *resid, double *direc, SistLinear_t *SL)
{
    double residTxresid = multiplicaVetores(resid, resid, SL->n);
    double alpha = 0;
    double *vetDirecxA = (double *)malloc(sizeof(double) * SL->n);
    double direcTxAxDirec = 0.0;

    // Percorre a matriz SL->A e o vetor de direcoes
    for (int i = 0; i < SL->n; ++i)
    {
        double soma = 0.0;
        for (int j = 0; j < SL->n; ++j)
        {
            // calcula o iésimo elemento do vetDirecxA (esse vetor é chamado de 'z' no livro M.Cristina C. Cunha)
            soma = soma + SL->A[i][j] * direc[j];

            // Teste para ver se não foi gerado um NaN ou um número infinito
            if (isnan(soma) || isinf(soma))
            {
                fprintf(stderr, "Erro soma(calcAlpha): %g é NaN ou +/-Infinito\n", soma);
                exit(1);
            }
        }
        // O iésimo elemento do vetor resultante recebe soma.
        vetDirecxA[i] = soma;
    }

    // Tendo o vetor resultante, agora é só multiplicar D^t * vetor_resultante.
    // O cálculo entre um vetor transposto e um vetor normal gera uma matriz de 1x1, ou seja, um número real.
    for (int i = 0; i < SL->n; ++i)
    {
        direcTxAxDirec += direc[i] * vetDirecxA[i];
        // Teste para ver se não foi gerado um NaN ou um número infinito
        if (isnan(direcTxAxDirec) || isinf(direcTxAxDirec))
        {
            fprintf(stderr, "Erro direcTxAxDirec(calcularDenominadorEscalarA): %g é NaN ou +/-Infinito\n", direcTxAxDirec);
            exit(1);
        }
    }

    alpha = residTxresid / direcTxAxDirec;

    free(vetDirecxA);
    return alpha;
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
 * @brief - Calcula o residuo   TODOOOOOOOOOOOOOOOO
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
    // alpha * A * p<k>
    for (int i = 0; i < n; ++i)
    {
        double somaI = 0;
        for (int j = 0; j < n; j++)
        {
            somaI += A[i][j] * p[j];
        }

        residuo[i] = residuo[i] - alpha * somaI;
    }
}

/**
 * @brief
 * beta = r^T * r / rAnt^T * rAnt
 * @param resid
 * @param residAnt
 * @param n
 * @return double
 */
double calcBeta(double *resid, double *residAnt, unsigned int n)
{
    double beta = multiplicaVetores(resid, resid, n) / multiplicaVetores(residAnt, residAnt, n);
    return beta;
}

/**
 * @brief Calcula a proxima direcao de busca p<k>
 *
 * p<k> = z<k> + beta<k-1> * p<k-1>
 *
 * @param proxDir - vetor da proxima direcao de busca
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

///////////////////////////////////////////////////////////
/**
 * @brief Metodo resolvedor de sistemas lineares pelo metodo de gradiente conjugado parece ok
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
    printf("AQUIIIII ");
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

    tempoPreCond = timestamp() - tempoResid;
    
    inicializaSol(x,SL->n);
    // residuo = b
    copiaVetor(SL->b, resid, SL->n);
    // direcao = residuo
    copiaVetor(resid, direc, SL->n);

    // for
    for (it = 0; it < maxIt; ++it)
    {
        double tIterInicio = timestamp();

        // CALCULA APLHA
        alpha = (multiplicaVetores(resid, resid, SL->n)) / (multiplicaVetores(vetorAux, direc, SL->n));
        // calcula novo x
        // x1 = x0 +alpha * p
        copiaVetor(x, xAnt, SL->n);
        calcX(x, xAnt, alpha, direc, SL->n);

        double normaMaxRel = normaMaxErroRelativo( x,xAnt, SL->n);
        fprintf(arqSaida, "# iter %d: %.15g \n", it, normaMaxRel);

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
        beta = calcBeta(resid, residAnt, SL->n);

        // calcula prox direc
        // P = R + b0 * pAnt
        copiaVetor(direc, direcAnt, SL->n);
        calcProxDirecBusca(direc, resid, beta, direcAnt, SL->n);

        tMedioIter += timestamp() - tIterInicio;
    }

    fprintf(arqSaida, "# residuo: || %.15g || \n", calcularNormaL2Residuo(resid, SL->n));

    // tempo final
    fprintf(arqSaida, "# Tempo PC: %.15g \n", tempoPreCond);
    tMedioIter = tMedioIter / it;
    fprintf(arqSaida, "# Tempo tMedioIter:: %.15g \n", tMedioIter);
    tempoResid = timestamp() - tempoResid;
    fprintf(arqSaida, "# Tempo residuo:: %.15g \n", tempoResid);
    free(resid);
    free(residAnt);
    free(direc);
    free(direcAnt);
    free(xAnt);
    free(x);

    return it;
}
