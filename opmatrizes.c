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

void transporMat(double **matA,double **matT, int n){
	for(int i=0;i < n ;++i)
		for(int j=0;j < n;++j){
			matT[i][j]=matA[j][i];  
		}
	return matT;
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