#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdlib>

const int MAXN = 505;
// typedef struct
// {
	// int row,col;
	// double *data; 
// }M;
class Vec
{
private: // nonstatic member
	int row;
	double dat[MAXN];
public: // static member
	Vec(int a);
	Vec();
	~Vec();
	friend Vec operator + (const Vec &A, const Vec &B); 
	friend Vec operator - (const Vec &A, const Vec &B);
	double &operator [] (const int x);
	double inprod(Vec &A,Vec &B);
};
Vec::Vec(int a): row(a)
{
}
Vec::~Vec()
{
}
Vec operator + (const Vec &A, const Vec &B)
{
	Vec ret(B.row);
	for (int i = 0;i < ret.row; i++)
	{
		ret.dat[i] = A.dat[i] + B.dat[i];
	}
	return ret;
}
Vec operator - (const Vec &A, const Vec &B)  // Why meiyou  Vec::
{
	Vec ret(B.row);
	for (int i = 0;i < ret.row; i++)
	{
		ret.dat[i] = A.dat[i] - B.dat[i];
	}
	return ret;
}
double &Vec::operator[] (const int x)
{
	return this->dat[x];
}
double Vec::inprod(Vec &A, Vec &B)
{
	double ret = 0;
	for (int i = 0;i < A.row;i++) 
		ret += A.dat[i] * B.dat[i];
	return ret;
}

class M
{
private:
	int row, col;
	double dat[MAXN][MAXN];

public:
	M(int a, int b): row(a), col(b) 
	{
	}
	M();
	~M();

	friend M operator * (const M &matrix_A, const M &matrix_B)
	{
		if (matrix_A.col != matrix_B.row)
			puts("Not compatible!"), exit(0);
		M ret(matrix_A.row,matrix_B.col);
		for(int i = 0;i < matrix_A.row;i++) 
			for(int j = 0;j < matrix_B.col;j++)
				for(int k = 0;k < matrix_A.col;k++)
				{
			 		ret.dat[i][j] += matrix_A.dat[i][k] * matrix_B.dat[k][j];
				}
		return ret;
	}
	friend M operator - (const M &matrix_A, const M &matrix_B)
	{
		if (matrix_A.row != matrix_B.row || matrix_A.col != matrix_B.col)
			puts("Not compatible!"), exit(0);
		M ret(matrix_A.row,matrix_A.col);
		for (int i = 0;i < matrix_A.row;i++)
			for (int j = 0;j < matrix_A.col;j++)
			{
				ret.dat[i][j] = matrix_A.dat[i][j] - matrix_B.dat[i][j];
			}
		return ret;
	}
	friend M operator + (const M &matrix_A, const M &matrix_B)
	{
		if (matrix_A.col != matrix_B.col || matrix_A.row!=matrix_B.row)
			puts("Not compatible!"), exit(0);
		M ret(matrix_A.row, matrix_A.col);
		for (int i=0;i<matrix_A.row;i++)
			for (int j=0;j<matrix_A.col;j++)
			{
				ret.dat[i][j] = matrix_A.dat[i][j] + matrix_B.dat[i][j];
			}
		return ret;
	}
	M operator ^ (long long exp)
	{
		if (this->col != this->row)
			puts("Not a square!"), exit(0);
		M Base = *this, ret = *this;
		for (exp--;exp;exp >>= 1ll, Base = Base * Base) 
			if (exp&1ll) 
				ret = ret * Base;
		return ret;
	}
	M Trans(M &matrix)
	{
		M ret(matrix.col,matrix.row);     
		for (int i=0;i<matrix.row;i++)
			for (int j=0;j<matrix.col;j++)
			{
				ret.dat[j][i] = matrix.dat[i][j];
			}
		return ret;
	}
	M Inv();
	
	// double inprod(const M &vecA, const M &vecB) { // inner product
	// 	double ret=0;
	// 	if (vecA.row != vecB.row) {puts("ERROR"); return 0;}
	// 	for (int i=0;i<vecA.row;i++)
	// 	{ 
	// 		ret += vecA.dat[i][0] * vecB.dat[i][0];
	// 	}
	// 	return ret;
	// }
	Vec proj(const Vec &vecA,const Vec &vecB, int n){
		double k = inprod(vecA, vecB) / inprod(vecB, vecB);
		Vec ret(n);
		for (int i = 0;i < n;i++)
		{
			ret[i] = k * vecB[i];
		}
		return ret;
	}
	void Schmidt (M &Q, M &R, const M matrix) {
		int n = matrix.row;
		Vec a[n], b[n]; // column vector
		for (int j = 0;j < n;j++)
			for (int i = 0;i < n;i++) 
			{
				a[j][i] = b[j][i] = matrix.dat[i][j];
			} //Partition matrix into n columns.   
		for (int i=0;i<n;i++)
			for (int j=0;j<i;j++) 
			{
				b[i] = b[i] - proj(a[i],b[j],n);
			}
		for (int j=0;j<n;j++)
		{
			double len = sqrt(Vec::inprod(b[j], b[j]));
			for (int i=0;i<n;i++)
			{
				Q.dat[i][j] = b[j].dat[i][0] / len;
			}
		} // merge n columns into Q
		// M* inv = InitMatrix(n,n); CopyMatrix(Q,inv);
		TransMatrix(inv);
		R = Trans(Q) * matrix;
		// MulMatrix(R,inv,matrix); // R= Q^-1 * mat = Q^T * mat
		// FreeMatrix(inv);
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


	M* InverseMatrix(M* Source);
	double detMatrix(M* Source); 
	void InitExpandMatrix(M* Source,M* expand);
	void getNewMatrix(M* res,M *Source);
};









M* InitMatrix(int row,int col);
M* _I(int n); //Identity Matrix
void Input(M* in,int row,int col); 
void Output();
int SizeMatrix(M* matrix);

// void MulMatrix(M* ret, M *matrix_A,M *matrix_B); // done 
// void DecMatrix(M* ret, M* matrix_A,M* matrix_B);
// void AddMatrix(M* ret, M *matrix_A,M *matrix_B);
// void PowerMatrix(M* ret, M* martix,long long  exp); 
// void TransMatrix(M* matrix);

// void Schmidt(M* Q,M* R, M* matrix);
// void Eig(M* matrix);
// double tr(M* matrix);
// void proj(M* ret, M* vecA, M* vecB );
// double inprod(M* vecA, M* vecB);

//���� 
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

}

int SizeMatrix(M *matrix)
{
	return matrix->row*matrix->col;
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

double tr(M* matrix){
	if (matrix->col != matrix->row) {puts("ERROR"); return 0;}
	double Ret=0;
	for (int i=0;i<matrix->col;i++)
		Ret+=Mat(matrix,i,i);//matrix->data[i*matrix->col+i];
	return Ret;
}

/* QR decomposition */
// double inprod(M* vecA, M* vecB) { // inner product
// 	double ret=0;
// 	if (vecA->row != vecB->row) {puts("ERROR"); return 0;}
// 	for (int i=0;i<vecA->row;i++)
// 	{ 
// 		ret+=Mat(vecA,i,0) * Mat(vecB,i,0);
// 	}
// 	return ret;
// }

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



// double det(M Source){//Сbug 2 1 0 1 
// 	int k=1,f=0;
// 	double t,ans=1.0,times;
// 	if( Source.col != Source.row) return pust("Not a square!"), 0;
// 		for(int i = 0, j = 0;i < Source.row-1 && j<Source.row - 1;i++, j++){
// 			if( Mat(Source,i,j) == 0 ){
// 				for(f = i; Mat(Source, f, j) == 0 && f<Source->row;++f);
// 				if(f == Source->row) return 0;
// 				else{
// 					for(int s=0;s < Source.col;++s){
// 						t=Mat(Source,i,s);
// 						Mat(Source,i,s) = Mat(Source,f,s);
// 						Mat(Source,f,s) = t;     
// 					}
// 					k=-k;
// 				}
// 			}
// 			printf("====== adjust ======\n");
// 			Output(Source);
// 			t=Mat(Source,i,j);
// 			for(int s=i+1;s<Source->row;++s){
// 				if(Mat(Source,s,j) != 0) 
// 				{
// 				times = Mat(Source,s,j)/t;
// 				for(int x=j;x<Source->col;++x)	
// 					Mat(Source,s,x) = Mat(Source,s,x) - Mat(Source,i,x) * times;
// 				}
// 			}	
// 		}
// 		printf("====== result ======\n");
// 			Output(Source);
// 	}
// 	for(int i=0,j=0; i<Source->row && j<Source->row; i++,j++){
// 		ans*=Mat(Source,i,j);
// 	return ans*k;
// }

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

