#include "utils.h"
#include "opmatrizes.h"

// performs dot product between two vectors
double produtoEscalar(double *vetA, double *vetB, int n)
{
	double somaTot=0;

	for (int i=0; i<n; ++i)
		somaTot+=vetA[i]*vetB[i];

	return somaTot;
}

// performs a*x+b operation with two vectors and a scalar
void aplicaX(double *dest, double a, double *x, double *b, int n)
{
	int i;
	for (i=0; i<n; ++i)
		dest[i]=a*x[i]+b[i];
}

/**
 * @brief Função para transpor uma matriz
 * 
 * @param matA - Matriz a ser transposta
 * @param matT - Matriz A^T
 * @param lin - numero de linhas
 * @param col - numero de colunas
 */
void transporMat(double **matA,double **matT, int lin, int col){
	for(int i=0;i < lin;++i)
		for(int j=0;j < col;++j){
			matT[j][i]=matA[i][j];  
		}
}

/**
 * @brief - Calcula o produto da matriz A com a matriz B
 * A * B = C
 * @param matA - Matriz A
 * @param matB - Matriz B
 * @param linA - Quantidade de linhas da matriz A
 * @param colA - Quantidade de colunas da matriz A
 * @param linB - Quantidade de linhas da matriz B
 * @param colB - Quantidade de colunas da matriz B
 * @return double** 
 */
double ** multMat(double ** matA, double ** matB, int linA, int colA,int linB, int colB)
{
	double **matrizResult = alocarMatriz(linA,colB); 
    for (int i = 0; i < linA; i++) {
        for (int j = 0; j < colB; j++) {
            matrizResult[i][j] = 0;

            for (int k = 0; k < linB; k++) {
                matrizResult[i][j] += matA[i][k] * matB[k][j];
            }
        } 
    }

	return matrizResult; 
}

double somarMatrizesTranspostas(double **a, double **b, int lin, int col)
{
	double c=0;
	for (int i = 0; i < lin; ++i)
		for (int j = 0; j < col; ++j)
			c+= a[i][j];
	
	for (int i = 0; i < col; ++i)
		for (int j = 0; j < lin; ++j)
			c+= b[i][j];
	
	return c; 
}

void geraMatrizIdentidade(double **m,int n){
	for(int i = 0; i < n;i++) {
		for(int j = 0; j < n;j++) {
			if(i == j) {
				m[i][j]=1;
			} else {
				m[i][j]=0;
			}
		}
	}
}