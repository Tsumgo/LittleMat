#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define Mat(mat,i,j) (mat)->data[(i)*(mat)->col+(j)]
#define FreeMatrix(matrix) { free(matrix->data); matrix->data = NULL;/* printf("释放成功\n"); */}

typedef struct
{
	int row,col;
	double *data; 
}M;
M* InitMatrix(int row,int col);
M* _I(int n); //Identity Matrix
void Input(); 
void Output(); 
void MulMatrix(M* ret, M *matrix_A,M *matrix_B); // done 
void SubtractMat(M* ret, M* matrix_A,M* matrix_B);
void AddMatrix(M* ret, M *matrix_A,M *matrix_B);
void PowerMatrix(M* ret, M* martix,long long  exp); 
void TransMatrix(M* matrix);
void Schmidt(M* Q,M* R, M* matrix);
void Eig(M* matrix);
int SizeMatrix(M* matrix);
double tr(M* matrix);


int n; long long k;
int main(){
freopen("in","r",stdin);

	scanf("%d%lld",&n,&k);
	M *A=InitMatrix(n,n), *Q=InitMatrix(n,n), *R=InitMatrix(n,n);
	Input(A,n,n);
	// Eig(A);
	// MulMatrix(A,A,A);
	Eig(A);
	// Schmidt(Q,R,A); // A = Q * R
	// Output(Q), Output(R);
	// MulMatrix(ans,A,A);
	// printf("A * A = \n"); 
	// Output(ans);
	// PowerMatrix(ans,A,k);
	// printf("A ^ k = \n");
	//  Output(ans);
	// Schmidt(ans,A);
	// Output(ans);
	// Output(_I(3));
	// Output(PowerMatrix(mat,k));  
}

int SizeMatrix(M *matrix)
{
	return matrix->row*matrix->col;
}

void CopyMatrix(M *matrix_A, M *matrix_B) //A --> B
{
	matrix_B->row = matrix_A->row;
	matrix_B->col = matrix_A->col;
	memcpy(matrix_B->data, matrix_A->data, SizeMatrix(matrix_A)*sizeof(double));
	///FreeMatrix(matrix_A);  // modified
}

void Input(M* matrix,int row,int col){
	for (int i=0;i<row;i++)
		for (int j=0;j<col;j++) 
			scanf("%lf",&Mat(matrix,i,j));
}
void Output(M* matrix){
	for (int i=0;i<matrix->row;i++){
		for (int j=0;j<matrix->col;j++)
			printf("%.5lf ",Mat(matrix,i,j));
		puts("");
	}
}

M* InitMatrix(int row,int col) {
	if (row<=0 || col <=0) return NULL;
	M* matrix = (M*)calloc(1,sizeof(M));
	matrix->col = col, matrix->row = row;
	matrix->data = (double*)calloc(row*col,sizeof(double)); // modified
	return matrix;
}
M* _I(int n){  // Identity Matrix
	M* matrix=InitMatrix(n,n); 
	for (int i=0;i<n;i++) Mat(matrix,i,i)=1;
	return matrix; 
}
double tr(M* matrix){
	if (matrix->col != matrix->row) {puts("ERROR"); return 0;}
	double Ret=0;
	for (int i=0;i<matrix->col;i++)
		Ret+=Mat(matrix,i,i);//matrix->data[i*matrix->col+i];
	return Ret;
}
/* + - * / */
void MulMatrix(M* ret, M *matrix_A,M *matrix_B){ 
	if (matrix_A->col != matrix_B->row) 
	{
		printf("ERROR\n");
		return;
	}
	M* ans = InitMatrix(matrix_A->row,matrix_B->col);
	for(int i=0;i<matrix_A->row;i++) 
		for(int j=0;j<matrix_B->col;j++)
		 	for(int k=0;k<matrix_A->col;k++)
			{
		 		Mat(ans,i,j) += Mat(matrix_A,i,k) * Mat(matrix_B,k,j);
			}
	CopyMatrix(ans, ret); FreeMatrix(ans);// ans --> ret
} 
void PowerMatrix(M* ret, M *Base,long long  exp){
	if (Base->col != Base->row)
	{
		printf("ERROR\n");
		return;
	}
	M* ans = _I(Base->col); 
	for (;exp;exp>>=1ll, MulMatrix(Base, Base, Base)) 
		if (exp&1)
		{
			MulMatrix(ans,ans,Base);
		}
	CopyMatrix(ans, ret); FreeMatrix(ans);
}
void DecMatrix(M* ret, M* matrix_A,M* matrix_B) {
	if (matrix_A->col != matrix_B->col || matrix_A->row!=matrix_B->row)
	{
		ret->col=ret->row=-1;
		return;
	}
	M* ans=InitMatrix(matrix_A->row,matrix_A->col);
	for (int i=0;i<matrix_A->row;i++)
	{
		for (int j=0;j<matrix_A->col;j++)
		{
			Mat(ans,i,j)=Mat(matrix_A,i,j)-Mat(matrix_B,i,j);
		}
	}
	CopyMatrix(ans, ret);  FreeMatrix(ans);
}
void AddMatrix(M* ret, M *matrix_A,M *matrix_B)
{
	if (matrix_A->row == matrix_B->row && matrix_A->col == matrix_B->col)
	{
		M *matrix_C = InitMatrix(matrix_A->row,matrix_A->col);
		for (int i=0;i<matrix_A->col;i++)
		{
			for (int j=0;j<matrix_A->row;j++)
			{
				matrix_C->data[i*matrix_C->row + j] = \
				matrix_A->data[i*matrix_A->row + j] + matrix_B->data[i*matrix_A->row + j];
			}
		}
		CopyMatrix(matrix_C, ret); FreeMatrix(matrix_C);
		// return matrix_C;   Modified
	}
	else 
	{
		printf("ERROR\n");
		// return NULL;
	}
}
void TransMatrix(M *matrix)			
{
	if (matrix->row == matrix->col)
	{
		M *matrixTemp = InitMatrix(matrix->row,matrix->col);     
		CopyMatrix(matrix,matrixTemp);

		for (int i=0;i<matrix->row;i++)
		{
			for (int j=0;j<matrix->col;j++)
			{
				matrix->data[i*matrix->row + j] = matrixTemp->data[j*matrix->row + i];
			}
		}
		FreeMatrix(matrixTemp);
	}
	else
	{
		printf("ERROR\n");
	}
}


/* QR decomposition */
double inprod(M* vecA, M* vecB) { // inner product
	double ret=0;
	if (vecA->row != vecB->row) {puts("ERROR"); return 0;}
	for (int i=0;i<vecA->row;i++)
	{ 
		ret+=Mat(vecA,i,0) * Mat(vecB,i,0);
	}
	return ret;
}
void proj(M* ret, M* vecA, M* vecB ){
	double k =inprod(vecA,vecB) 
			/ inprod(vecB, vecB);
	M* ans = InitMatrix(n,1);
	for (int i=0;i<vecA->row;i++)
	{
		Mat(ans,i,0)=k*Mat(vecB,i,0);
	}
	CopyMatrix(ans,ret);  FreeMatrix(ans);
}
void Schmidt (M* Q, M* R, M* matrix) {
	int n=matrix->row;
	M *a[n], *b[n]; 
	for (int j=0;j<n;j++)
	{
		a[j] = InitMatrix(n, 1), b[j] = InitMatrix(n,1);
		for (int i=0;i<n;i++) Mat(b[j],i,0) = Mat(a[j],i,0) = Mat(matrix,i,j);
	} //Partition matrix into n columns.
	M* temp = InitMatrix(n,1);
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<i;j++) 
		{
			proj(temp,a[i],b[j]); // temp= ||a[i],b[j]|| / ||b[j], b[j]|| * b[j]
			DecMatrix(b[i], b[i], temp); // b[i] = b[i] - temp;
		}
	}  
	FreeMatrix(temp);
	for (int j=0;j<n;j++)
	{
		double len = sqrt(inprod(b[j],b[j]));
		for (int i=0;i<n;i++)
		{
			Mat(Q,i,j)=Mat(b[j],i,0)/len;
		}
		FreeMatrix(a[j]); FreeMatrix(b[j]);
	} // merge n columns into Q
	M* inv = InitMatrix(n,n); CopyMatrix(Q,inv);
	TransMatrix(inv);
	MulMatrix(R,inv,matrix); // R= Q^-1 * mat = Q^T * mat
	FreeMatrix(inv);
}
void Eig(M* matrix)
{
	int n=matrix->col;
	M* Q=InitMatrix(n,n), *R=InitMatrix(n,n);
	M* inv=InitMatrix(n,n);
	for (int i=1;i<=1000;i++)
	{
		Schmidt(Q,R,matrix);
		MulMatrix(matrix,R,Q);
	}
	for (int i=0;i<n;i++)
		for (int j=i+1;j<n;j++) Mat(R,i,j)=0;
	Output(R);
}



