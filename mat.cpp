#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdlib>
#define abs(x) ((x) > 0 ? (x): -(x))


const int MAXN = 500;

void swap(double &x,double &y) 
{
	double tmp = x;
	x=y, y=tmp;
	// x^=y^=x^=y;
}


class Vec
{
protected: // nonstatic member
	double *dat;
public: // static member
	int row;
	Vec(int a):row(a)
	{
		dat = (double*) calloc(a, sizeof (double));
	}
	Vec(){}
	~Vec()
	{
		
	}
	void input()
	{
		for (int i = 0;i < row;i++) 
			scanf("%lf", &dat[i]);
	}
	void print()
	{
		for (int i = 0;i < row;i++) 
			printf("%.2lf ",dat[i]);
		puts("");
	}
	double &operator [] (const int x) // operator [] 必须是成员函数？
	{
		return this->dat[x];
	}
};
Vec operator + (Vec &A, Vec &B) // 全局函数
{
	Vec ret(B.row);
	for (int i = 0;i < ret.row; i++)
	{
		ret[i] = A[i] + B[i];
	}
	return ret;
}
Vec operator - (Vec &A, Vec B)
{
	Vec ret(B.row);
	for (int i =0 ;i< ret.row;i++)
	{
		ret[i] = A[i] - B[i];
	}
	return ret;
}
double inprod(Vec &A, Vec &B)
{
	double ret = 0;
	for (int i = 0;i < A.row;i++) 
		ret += A[i] * B[i];
	return ret;
}
Vec proj(Vec &vecA, Vec &vecB)
{
	double k = inprod(vecA, vecB) 
			/ inprod(vecB, vecB);
	int n=vecA.row;
	Vec ret(n); 
	for (int i = 0;i < n;i++)
	{
		ret[i] = k * vecB[i];
	}
	return ret;
}
class M
{
protected:
	double dat[MAXN][MAXN];
public:
	int row, col;
	M(int a, int b)
	{
		row = a, col = b;
	}
	M(){}
	~M(){}
	M Trans();
	friend M operator * (const M &matrix_A, const M &matrix_B);
	friend M operator - (const M &matrix_A, const M &matrix_B);
	friend M operator + (const M &matrix_A, const M &matrix_B);
	friend void Schmidt (M &Q, M &R, const M matrix); 
	void Eig(M matrix)
	{
		int n=matrix.col;
		M Q(n,n), R(n,n);
		// if(detMatrix(matrix) == 0 )
		// {
		// 	printf("不好分解");
		// 	return; 
		// }
		for (int i=1;i<=1000;i++)
		{
			Schmidt(Q,R,matrix);
			matrix = R * Q;
		}
		for (int i=0;i<n;i++)
			for (int j=i+1;j<n;j++) R.dat[i][j]=0;
		R.Output();
	}

	double tr();
	M Elimination(M Source);
	double Det();
	M Inv();
	void Output();
	void Input(int row, int col);
};

void M::Input(int row,int col){
	for (int i=0;i<row;i++)
		for (int j=0;j<col;j++) 
			scanf("%lf",&this->dat[i][j]);
}
void M::Output(){
	for (int i=0;i<this->row;i++){
		for (int j=0;j<this->col;j++)
			printf("%5.1lf ",this->dat[i][j]);
		printf("\n");
	}
}
M operator * (const M &matrix_A, const M &matrix_B)
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
M operator - (const M &matrix_A, const M &matrix_B)
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
M operator + (const M &matrix_A, const M &matrix_B)
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
M operator ^ (const M &matrix, long long exp)
{
	if (matrix.col != matrix.row)
		puts("Not a square!"), exit(0);
	M Base = matrix, ret = matrix;
	for (exp--;exp;exp >>= 1ll, Base = Base * Base) 
		if (exp&1ll) 
			ret = ret * Base;
	return ret;
}
M M::Trans()
{
	M ret(col, row);     
	for (int i=0;i<row;i++)
		for (int j=0;j<col;j++)
		{
			ret.dat[j][i] = dat[i][j];
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
			b[i] = b[i] - proj(a[i],b[j]);
		}
	for (int j = 0;j < n;j++)
	{
		double len = sqrt(inprod(b[j], b[j]));
		for (int i=0;i<n;i++)
		{
			Q.dat[i][j] = b[j][i] / len;
		}
	} // merge n columns into Q
	R = Q.Trans() * matrix;
}
double M::tr(){
	if (col != row) {puts("Not a square!"); return 0;}
	double ret = 0;
	for (int i = 0;i < col;i++)
		ret += dat[i][i];
	return ret;
}
double M::Det()
{  // 这里会修改Source的值，就不引用了。 
	int k = 1, f = 0;
	double t, ans = 1.0, times;
	if( col != row) return puts("Not a square!"), 0;
	M tmp = Elimination(*this);
		// for(int i = 0, j = 0;i < row-1 && j < row - 1;i++, j++){
		// 	if( dat[i][j] == 0 ){
		// 		for(f = i; dat[f][j] == 0 && f<row;++f);
		// 		if(f == row) return 0;
		// 		else{
		// 			for(int s=0;s < col;++s){
		// 				t=dat[i][s]; //Mat(Source,i,s);
		// 				dat[i][s] = dat[f][s]; //Mat(Source,f,s);
		// 				dat[f][s] = t;
		// 			}
		// 			k=-k;
		// 		}
		// 	}
		// 	printf("====== adjust ======\n");
		// 	Output();
		// 	t = dat[i][j];
		// 	for(int s = i + 1;s < row;++s){
		// 		if(dat[s][j] != 0) 
		// 		{
		// 		times = dat[s][j]/t;
		// 		for(int x = j;x < col;++x)	
		// 			dat[s][x] = dat[s][x] - Source.dat[i][x] * times;
		// 		}
		// 	}	
		// }
	for (int i = 0; i < row;i++) ans*=dat[i][i];
	return ans;
}
M M::Elimination(M expand)
{ // 阶梯型 阶梯头为1。
		int f, i = 0, j = 0;
		while (i < expand.row && j < expand.col)
		{
			f = i;
			for (int k = 0;k < expand.row; k++)
				if (abs(expand.dat[k][j]) > abs(expand.dat[f][j]))
				{
					f = k;
				}
			for (int k = 0;k < expand.col;k++)
			{
				swap(expand.dat[i][k], expand.dat[f][k]);
			}
			if (expand.dat[i][j] == 0)
			{
				j++;
				continue;
			}
			for (int k = expand.col - 1;k >= j;k--)
			{ 
				expand.dat[i][k] /= expand.dat[i][j];
			}
			for (int k = i+1; k < expand.row; ++k)
				for (int l = expand.col - 1;l >= j;l--)
				{
					expand.dat[k][l] -= expand.dat[k][j];
				}
			i++;
		}
		return expand;
}
M  M::Inv()
{
	if (col != row) 
		puts("Peculiar!"), exit(0);
	int n = col;
	M ans(n, n);
	M expand(n, n*2);
	for (int i = 0;i < n;i++)
		for  (int j = 0;j < n*2;j++)
		{
			if (j < n) expand.dat[i][j] = dat[i][j];
			else expand.dat[i][j] = (j-n == i);
		}
	expand = Elimination(expand);
	for (int i = 0; i < n;i++) 
		if (expand.dat[i][i] == 0) 
			puts("Peculiar!"), exit(0);
	for (int i = 0;i < n;i++)
		for (int j = n;j < n*2;j++)
			ans.dat[i][j-n] = expand.dat[i][j];
	return ans;
}
 
int n,m; long long k;
int shit()
{
	// M A(1,1);
	// printf("%d\n",sizeof (A));

}
int main(void) 
{
	freopen("in","r",stdin);
	scanf("%d %d",&n,&m);
	shit();
	// M A(n,m), B(n,m);
	// M B(n,n);
	// A.input(), B.input();
	// A.print();
	// B.print();
	// (A+B).print();
}
// int main( int argc,char *argv[]){
/*
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
		M equations(p,q+1);
		printf("====== init_equations =====\n"); //���뷽���� 
		equations.Input(equations.row,equations.col);
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
*/
// }