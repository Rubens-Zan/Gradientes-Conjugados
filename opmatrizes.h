double produtoEscalar(double *vetA, double *vetB, int n);
void transporMat(double **matA,double **matT, int lin, int col);
void geraMatrizIdentidade(double **m,int n); 
void multMat(double ** matA, double ** matB, int linA, int colA,int linB, int colB,double **matrizResult);
double formMultMatrizesGeraValor(double ** matA, double ** matB, int linA, int colA,int linB, int colB);


double multiplicaVetor_Vetor (double *a, double *b, unsigned int N);
void multiplicaInteiro_Vetor (double mult, double *v, double *vetorAux, unsigned int N);
void somaVetor (double *a, double *b, double *vetorFinal, unsigned int N);
void copiaVetor (double *a, double *b, unsigned int N);
void multiplicaMatriz_Vetor(double **matriz, double *v, double *z, unsigned int N);
void subtraiVetor (double *a, double *b, double *vetorFinal, unsigned int N);