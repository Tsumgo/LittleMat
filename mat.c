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

//�������� 
M* InitMatrix(int row,int col);
M* _I(int n); //Identity Matrix
void Input(M* in,int row,int col); 
void Output();
int SizeMatrix(M* matrix);

//�������� 
void MulMatrix(M* ret, M *matrix_A,M *matrix_B); // done 
void DecMatrix(M* ret, M* matrix_A,M* matrix_B);
void AddMatrix(M* ret, M *matrix_A,M *matrix_B);
void PowerMatrix(M* ret, M* martix,long long  exp); 
void TransMatrix(M* matrix);

//��������������ֵ 
void Schmidt(M* Q,M* R, M* matrix);
void Eig(M* matrix);
double tr(M* matrix);
void proj(M* ret, M* vecA, M* vecB );
double inprod(M* vecA, M* vecB);

//���� 
M* InverseMatrix(M* Source);
double detMatrix(M* Source); 
void InitExpandMatrix(M* Source,M* expand);
void getNewMatrix(M* res,M *Source);

//��˹��Ԫ�ⷽ���� 
M* elimination(M* Source);
 
int n; long long k;
int main( int argc,char *argv[]){
 	//freopen("in","r",stdin); //���ļ�����
	
	if(*argv[1] != '-'  ) //Сbug���� -- 
	{
		printf("ע�������ʽ��%s -sth.",argv[0]);
		return 0;
	}
	if( *(argv[1]+1) == 'h')
	{
	printf("====== help =======\n");
	printf("-+: ����ӷ�\n");
	printf("--: �������\n");
	printf("-*: ����˷�\n");
	printf("-m����������\n");
	printf("-t: ����ת��\n");
	printf("-d��������ʽ\n");
	printf("-i: �������\n"); 	
	printf("-s: ������ֵ\n");
	printf("-x: �ⷽ����\n");
	return 0;
	}
	if(*(argv[1]+1) == 'x')
	{
		int p,q;
		printf("OVO�����м������̣�");
		scanf("%d",&p);
		printf("OVO�����м���δ֪����");
		scanf("%d",&q);
		M* equations = InitMatrix(p,q+1);
		printf("====== init_equations =====\n"); //���뷽���� 
		Input(equations,equations->row,equations->col);
		elimination(equations);
		for(int i=0;i<equations->row;++i)
			printf("x%d = %.1lf, ",i+1,Mat(equations,i,equations->col-1));
		return 0;	
	}
	printf("�������Ľ�����");
	scanf("%d",&n); 
	M *A=InitMatrix(n,n), *Q=InitMatrix(n,n), *R=InitMatrix(n,n);
	if(*(argv[1]+1) == 's') //��������������ֵ 
	{
		printf("====== init_matrix ======\n");
		Input(A,n,n);
		printf("====== result_matrix ======\n");
		Eig(A);
	}
	if( *(argv[1]+1) == '+')
	 {
 			 printf("====== first_matrix ======\n");
 			 Input(Q,n,n);
 			 printf("====== second_matrix ======\n");
 			 Input(R,n,n);
 			 AddMatrix(A,Q,R);
 			 printf("====== result_matrix ======\n");
 			 Output(A);
	 }
	 if( *(argv[1]+1) == '-')
	 {
	 		 printf("====== first_matrix ======\n");
 			 Input(Q,n,n);
 			 printf("====== second_matrix ======\n");
 			 Input(R,n,n);
 			 DecMatrix(A,Q,R);
 			 printf("====== result_matrix ======\n");
 			 Output(A);
	 }
	 if( *(argv[1]+1) == '*')
	 {
	 		 printf("====== first_matrix ======\n");
 			 Input(Q,n,n);
 			 printf("====== second_matrix ======\n");
 			 Input(R,n,n);
 			 MulMatrix(A,Q,R);
 			 printf("====== result_matrix ======\n");
 			 Output(A);
	 }
	 if(*(argv[1]+1) == 'd')
	 {
	    	printf("====== init_matrix ======\n");	
	 		Input(A,n,n);
	 		printf("the determinant is %lf",detMatrix(A));
	 } 
	if(*(argv[1]+1) == 'i')
	{
		 printf("====== init_matrix ======\n");
		 Input(A,n,n);

		 Output(InverseMatrix(A));
	}
	if(*(argv[1]+1) == 'p')
	{
		long long k;
		printf("====== init_matrix ======\n");	
		Input(A,n,n);
		printf("====== mutiple_times ======\n");
		scanf("%lld",&k);
		PowerMatrix(Q,A,k);
		printf("====== result_matrix =====\n");
		Output(Q);
	}
	if(*(argv[1]+1) == 't')
	{
		printf("====== init_matrix ======\n");
		Input(A,n,n);
		TransMatrix(A);
		printf("====== result_matrix =====\n");
		Output(A);
	}

	//Input(A,n,n);
	//Eig(A);
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
			printf("%5.1lf ",Mat(matrix,i,j));
		printf("\n");
	}
}

M* InitMatrix(int row,int col) {
	if (row<=0 || col <=0) return NULL;						
	M* matrix = (M*)calloc(1,sizeof(M));                    //calloc��malloc�������� �� ��calloc�Կռ���г�ʼ����ȫΪ0 
	matrix->col = col, matrix->row = row;
	matrix->data = (double*)calloc(row*col,sizeof(double)); // modified
	return matrix;											//����ֱ�Ϊ�ṹ��ָ��ͽṹ���е�dataָ�뿪�ռ䣬��ֹ��ʧ 
}


M* _I(int n){  // ����һ��n�׵ĵ�λ���� 
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
	
	CopyMatrix(ans, ret); FreeMatrix(ans);// ans --> ret ������ʵ����һ��С�Ż�������ans���浽һ��ָ���ĵط� ret���ͼ�ʱ�ͷ� 
										 // ��������������������ʵ�֣���ֹMulMatrix����һֱ���ɾ���ռ�ݿռ� 
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
					Mat(matrix_C,i,j)=Mat(matrix_A,i,j)+Mat(matrix_B,i,j);
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
	if(detMatrix(matrix) == 0 )
	{
		printf("�þ���Ϊ�������ʱ�޷������");
		return; 
	}
	for (int i=1;i<=1000;i++)
	{
		Schmidt(Q,R,matrix);
		MulMatrix(matrix,R,Q);
	}
	for (int i=0;i<n;i++)
		for (int j=i+1;j<n;j++) Mat(R,i,j)=0;
	Output(R);
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



double detMatrix(M *Source){//Сbug 2 1 0 1 
	int k=1,f=0;
	double t,ans=1.0,times;
	if( Source->col != Source->row) return 0;
	else{
		for(int i=0,j=0;i<Source->row-1 && j<Source->row-1;i++,j++){
			if( Mat(Source,i,j) == 0 ){
				for(f=i;Mat(Source,f,j) == 0 && f<Source->row;++f);
				if(f == Source->row) return 0;
				else{
					for(int s=0;s<Source->col;++s){
						t=Mat(Source,i,s);
						Mat(Source,i,s) = Mat(Source,f,s);
						Mat(Source,f,s) = t;     
					}
					k=-k;
				}
			}
			printf("====== adjust ======\n");
			Output(Source);
			t=Mat(Source,i,j);
			for(int s=i+1;s<Source->row;++s){
				if(Mat(Source,s,j) != 0) 
				{
				times = Mat(Source,s,j)/t;
				for(int x=j;x<Source->col;++x)	
					Mat(Source,s,x) = Mat(Source,s,x) - Mat(Source,i,x) * times;
				}
			}	
		}
		printf("====== result ======\n");
			Output(Source);
	}
	for(int i=0,j=0; i<Source->row && j<Source->row; i++,j++){
		ans*=Mat(Source,i,j);
	}
	return ans*k;
}

M*  InverseMatrix(M *Source){
	M* ans=InitMatrix(Source->row,Source->col);
	int f;
	double t,times;
		{
		M* expand = InitMatrix (Source->row,Source->col*2);
		InitExpandMatrix(Source,expand);
		for(int i=0,j=0;i<Source->row && j<Source->col;i++,j++)
		{
			if( Mat(expand,i,j) == 0 && i<Source->row-1 && j<Source->col-1){
				for(f=i+1;Mat(expand,f,j) == 0 && f<Source->row;++f);
				{
					for(int s=0;s<expand->col;++s)
					{
						t=Mat(expand,i,s);
						Mat(expand,i,s) = Mat(expand,f,s);
						Mat(expand,f,s) = t;     
					}
					if(f == Source->row); 
					{
					printf("�þ��󲻿��棡"); 
					return NULL; 
					}
				}	
			if( Mat(expand,i,j) == 0 && i==Source->row-1 && j==Source->col-1)
				{
				printf("�þ��󲻿��棡");
				return NULL;
				}
			}
			t=Mat(expand,i,j);
			for(int s = 0; s<expand->row; ++s)
			{
				if( s == i)
				{
					continue;
				}
				if( Mat(expand,s,j) != 0 ) 
				{
					times = Mat(expand,s,j)/t;
					for(int x=j;x<expand->col;++x)	
						Mat(expand,s,x) = Mat(expand,s,x) - Mat(expand,i,x) * times;
				}
			} 
		}
		for(int i=0,j=0;i<Source->row && j<expand->col;i++,j++)
		{
			t = Mat(expand,i,j);
			for(int x=0;x<expand->col;x++)
				Mat(expand,i,x) = Mat(expand,i,x)/t;
		}
		printf("======expand_matrix======\n");
		Output(expand);
	getNewMatrix(ans,expand); 
		printf("====== inverse_matrix =====\n");
	return ans;
}
}

void getNewMatrix(M* res,M *Source){
	for (int i=0;i<Source->row;++i){
		for (int j=res->col;j<Source->col;++j){
			Mat(res,i,j-res->col) = Mat(Source,i,j);
		}
}
}

void InitExpandMatrix(M* Source,M* expand){
	int i,j;
	for(i=0;i<expand->row;++i){
		for(j=0;j<expand->col;++j){
			if(j<Source->col) Mat(expand,i,j) = Mat(Source,i,j);
			else if( i == j-Source->col) Mat(expand,i,j) = 1; 
		}
	}
}

M* elimination(M* expand){
		double t,times;
		int f;
		for(int i=0,j=0;i<expand->row && j<expand->col;i++,j++)
		{
			if( Mat(expand,i,j) == 0 && i<expand->row-1 && j<expand->col-1)
			{
				for(f=i+1;Mat(expand,f,j) == 0 && f<expand->row;++f);
				{
					for(int s=0;s<expand->col;++s){
						t=Mat(expand,i,s);
						Mat(expand,i,s) = Mat(expand,f,s);
						Mat(expand,f,s) = t;     
					}
				}
				if(f == expand->row)
				{
					printf("����������������޽⣡"); 
					return NULL;
				}
			}
			if( Mat(expand,i,j) == 0 && i==expand->row-1 && j==expand->col-1)
				{
				printf("����������������޽⣡");
				return NULL;
				}
			t=Mat(expand,i,j);
			for(int s = 0; s<expand->row; ++s)
			{
				if( s == i)
				{
					continue;
				}
				if( Mat(expand,s,j) != 0 ) 
				{
					times = Mat(expand,s,j)/t;
					for(int x=j;x<expand->col;++x)	
						Mat(expand,s,x) = Mat(expand,s,x) - Mat(expand,i,x) * times;
				}
			} 
		}
		for(int i=0,j=0;i<expand->row && j<expand->col;i++,j++)
		{
			t = Mat(expand,i,j);
			for(int x=0;x<expand->col;x++)
				Mat(expand,i,x) = Mat(expand,i,x)/t;
		}
		Output(expand);
}

