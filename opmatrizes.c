#include "utils.h"
#include "opmatrizes.h"
#include "sislin.h"
#include "aritmetica.h"

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
 * @param matrizResult - Matriz resultante da multiplicação
 * @return double** 
 */
void multMat(double ** matA, double ** matB, int linA, int colA,int linB, int colB,double **matrizResult)
{
    for (int i = 0; i < linA; i++) {
        for (int j = 0; j < colB; j++) {
            matrizResult[i][j] = 0;

            for (int k = 0; k < linB; k++) {
                matrizResult[i][j] += matA[i][k] * matB[k][j];
            }
        } 
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
double formMultMatrizesGeraValor(double ** matA, double ** matB, int linA, int colA,int linB, int colB)
{
	double **matrizResult = alocarMatriz(linA+1,colB+1); 

	multMat(matA,matB,linA,colA, linB,colB,matrizResult);
	double valorFormatado = matrizResult[0][0];


	printf("VALOR FORMATADO :: %f\n", valorFormatado);
	prnMat(matrizResult,linA, colB); 
	
	liberarMatriz(matrizResult);
	return valorFormatado; 
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




//////////////////////////

//Multiplica Matriz com vetor
/**
 * @brief Multiplica uma matriz e um vetor
 * @param matriz(n,n) - Matriz a ser multiplicada
 * @param v - Vetor
 * @param z - Vetor que vai receber o produto
 * @param N - Dimensao da matriz
 */
void multiplicaMatriz_Vetor(double **matriz, double *v, double *z, unsigned int N) {
    double sum;
    
	for (int i = 0; i < N; i++) {
        double soma = 0;
        
		for (int j = 0; j < N; j++) {
            soma += matriz[i][j] * v[j];
        }
        
		z[i] = soma;
    }
}

//Multiplica vetor com vetor e retorna o resultado
/**
 * @brief - Multiplica dois vetores retornando a soma
 * R = a * b
 * @param a 
 * @param b 
 * @param N 
 * @return double 
 */
double multiplicaVetor_Vetor (double *a, double *b, unsigned int N) {
    double result = 0;
    int i;

    for (i = 0; i < N; i++) {
        result += a[i]*b[i];
    }
    
	return result;

}

//Muiltiplica um número pelo vetor
/**
 * @brief 
 * VETOR AUX = MULT * V
 * @param mult - MULTIPLO 
 * @param v 
 * @param vetorAux 
 * @param N 
 */
void multiplicaInteiro_Vetor (double mult, double *v, double *vetorAux, unsigned int N) {
	int i;

	for (i = 0; i < N; i++) {
		vetorAux[i] = mult*v[i];
	}
}

//Soma 2 vetores
/**
 * @brief - SOMA 2 VETORES E ARMAZENA EM VF
 * VF = A + B
 * @param a 
 * @param b 
 * @param vetorFinal 
 * @param N 
 */
void somaVetor (double *a, double *b, double *vetorFinal, unsigned int N) {
	int i;

	for (i = 0; i < N; i++) {
		vetorFinal[i] = a[i]+b[i];
	}
}

//Subtrai 2 vetores
/**
 * @brief - Subtrai dois vetores
 * VF = A - B
 * @param a 
 * @param b 
 * @param vetorFinal 
 * @param N 
 */
void subtraiVetor (double *a, double *b, double *vetorFinal, unsigned int N) {
	int i;

	for (i = 0; i < N; i++) {
		vetorFinal[i] = a[i]-b[i];
	}
}

//Copia um vetor
void copiaVetor (double *a, double *b, unsigned int N) {
	int i;
	
	for (i = 0; i < N; i++) {
		b[i] = a[i];
	}
}
