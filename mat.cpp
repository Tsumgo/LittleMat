#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstdlib>
#include "Matrix.hpp"
#define abs(x) ((x) > 0 ? (x) : -(x))

const int MODE = 1e9 + 7;

void swap(double &x, double &y)
{
	double tmp = x;
	x = y, y = tmp;
}

Vec::Vec(int a) : row(a) { dat = new double[a]; }
Vec::Vec() {}
Vec::~Vec() {}
void Vec::input()
{
	for (int i = 0; i < row; i++)
		scanf("%lf", &dat[i]);
}
void Vec::print()
{
	for (int i = 0; i < row; i++)
		printf("%.2lf ", dat[i]);
	puts("");
}
double &Vec::operator[](const int x) // operator [] 必须是成员函数？
{
	return this->dat[x];
}
Vec operator+(Vec &A, Vec &B) // 全局函数
{
	Vec ret(B.row);
	for (int i = 0; i < ret.row; i++)
	{
		ret[i] = A[i] + B[i];
	}
	return ret;
}
Vec operator-(Vec &A, Vec B)
{
	Vec ret(B.row);
	for (int i = 0; i < ret.row; i++)
	{
		ret[i] = A[i] - B[i];
	}
	return ret;
}
double inprod(Vec &A, Vec &B)
{
	double ret = 0;
	for (int i = 0; i < A.row; i++)
		ret += A[i] * B[i];
	return ret;
}
Vec proj(Vec &vecA, Vec &vecB)
{
	double k = inprod(vecA, vecB) / inprod(vecB, vecB);
	int n = vecA.row;
	Vec ret(n);
	for (int i = 0; i < n; i++)
	{
		ret[i] = k * vecB[i];
	}
	return ret;
}
/*Vector*/
/*Matrix*/
Matrix *Matrix::creatmatrix(int row, int col)
{
	Matrix *mat = new Matrix;
	mat->row = row, mat->col = col;
	mat->dat = new double *[row]; // Why new (double *) [row] is wrong?
	for (int i = 0; i < row; i++)
		mat->dat[i] = new double[col];
	return mat;
}
void Matrix::recycle(Matrix *matrix)
{
	for (int i = 0; i < matrix->row; i++)
		delete[] matrix->dat[i];
	delete[] matrix->dat;
	delete matrix;
}
void Matrix::Input(int row, int col)
{
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			scanf("%lf", &this->dat[i][j]);
}
void Matrix::print()
{
	for (int i = 0; i < this->row; i++)
	{
		for (int j = 0; j < this->col; j++)
			printf("%5.3lf ", this->dat[i][j]);
		printf("\n");
	}
	puts("---------");
}
void Matrix::multi(Matrix *C, const Matrix *matrix_A, const Matrix *matrix_B)
{
	Matrix *temp = creatmatrix(matrix_A->row, matrix_B->col), *matrix_B1 = creatmatrix(matrix_B->col, matrix_B->row);
	if (matrix_A->col != matrix_B->row)
		puts("Not compatible!"), exit(0);
	Trans(matrix_B1, matrix_B); // 转置后乘
	for (int i = 0; i < matrix_A->row; i++)
		for (int j = 0; j < matrix_B->col; j++)
		{
			double sum = 0;
			for (int k = 0; k < matrix_A->col; k++)
			{
				sum += matrix_A->dat[i][k] * matrix_B1->dat[j][k]; // modified
			}
			temp->dat[i][j] = sum;
		}
	for (int i = 0; i < matrix_A->row; i++)
		for (int j = 0; j < matrix_B->col; j++)
			C->dat[i][j] = temp->dat[i][j];
	recycle(temp), recycle(matrix_B1);
}
void Matrix::strassen(Matrix *C, const Matrix *A, const Matrix *B, int size)
{
	if (size <= 64)
	{
		multi(C, A, B);
		return;
	}

	// 分块
	//
}
void Matrix::sub(Matrix *C, const Matrix *matrix_A, const Matrix *matrix_B)
{
	if (matrix_A->row != matrix_B->row || matrix_A->col != matrix_B->col)
		puts("Not compatible!"), exit(0);
	for (int i = 0; i < matrix_A->row; i++)
		for (int j = 0; j < matrix_A->col; j++)
		{
			C->dat[i][j] = matrix_A->dat[i][j] - matrix_B->dat[i][j];
		}
}
void Matrix::add(Matrix *C, const Matrix *matrix_A, const Matrix *matrix_B)
{
	if (matrix_A->col != matrix_B->col || matrix_A->row != matrix_B->row)
		puts("Not compatible!"), exit(0);
	for (int i = 0; i < matrix_A->row; i++)
		for (int j = 0; j < matrix_A->col; j++)
		{
			C->dat[i][j] = matrix_A->dat[i][j] + matrix_B->dat[i][j];
		}
}
void Matrix::quickpow(Matrix *C, Matrix *Base, long long exp)
{ // Base 会被改变
	if (Base->col != Base->row)
		puts("Not a square!"), exit(0);
	for (int i = 0; i < Base->col; i++)
		C->dat[i][i] = 1;
	C->col = C->row = Base->col;
	for (; exp; exp >>= 1ll, multi(Base, Base, Base))
		if (exp & 1ll)
			multi(C, C, Base);
}
void Matrix::Trans(Matrix *result, const Matrix *matrix)
{
	for (int i = 0; i < matrix->row; i++)
		for (int j = 0; j < matrix->col; j++)
		{
			result->dat[j][i] = matrix->dat[i][j];
		}
}
void Matrix::Schmidt(Matrix *Q, Matrix *R, const Matrix *matrix)
{
	int n = matrix->row;
	Vec a[n], b[n]; // column vector
	// Matrix temp(matrix->row, matrix->col);
	for (int j = 0; j < n; j++)
		for (int i = 0; i < n; i++)
		{
			a[j].row = b[j].row = n;
			a[j][i] = b[j][i] = matrix->dat[i][j];
		} // Partition matrix into n columns->
	for (int i = 0; i < n; i++)
		for (int j = 0; j < i; j++)
		{
			b[i] = b[i] - proj(a[i], b[j]);
		}
	for (int j = 0; j < n; j++)
	{
		double len = sqrt(inprod(b[j], b[j]));
		for (int i = 0; i < n; i++)
		{
			Q->dat[i][j] = b[j][i] / len;
		}
	} // merge n columns into Q
	Matrix *Q_trans = creatmatrix(matrix->row, matrix->col);
	Trans(Q_trans, Q);
	multi(R, Q_trans, matrix);
	recycle(Q_trans);
}
void Matrix::Eig()
{
	Matrix *matrix = this;
	int n = col;
	Matrix *Q = creatmatrix(n, n), *R = creatmatrix(n, n);
	for (int i = 1; i <= 1000; i++)
	{
		Schmidt(Q, R, matrix);
		multi(matrix, R, Q);
	}
	for (int i = 0; i < n; i++)
		for (int j = i + 1; j < n; j++)
			R->dat[i][j] = 0;
	R->print();
}
double Matrix::trace(Matrix *source)
{
	if (source->col != source->row)
	{
		puts("Not a square!");
		return 0;
	}
	double ret = 0;
	for (int i = 0; i < source->col; i++)
		ret += source->dat[i][i];
	return ret;
}
double Matrix::Det(Matrix *source) //
{
	double ans;
	int i = 0, j = 0, f, cnt = 0;
	if (source->col != source->row)
		return puts("Not a square!"), 0;
	Matrix mat = *source; // 栈区， 会自行删除
	while (i < mat.row && j < mat.col)
	{
		f = i;
		for (int k = i + 1; k < mat.row; k++)
			if (abs(mat.dat[k][j]) > abs(mat.dat[f][j]))
				f = k;
		if (f != i)
			for (int k = 0; k < mat.col; k++)
			{
				swap(mat.dat[i][k], mat.dat[f][k]);
			}
		else
			cnt++;
		if (mat.dat[i][j] == 0)
			return 0;
		for (int k = i + 1; k < mat.row; ++k)
			mat.eli(i, j, k);
		i++, j++;
	}
	ans = cnt % 2 ? -1 : 1;
	for (int i = 0; i < mat.row; i++)
		ans *= mat.dat[i][i];
	return ans;
}
void Matrix::eli(int x, int y, int z) // C[z] -= C[x] * dat[z][y]; 用x行消去z行，使z,y处变成阶梯头。
{
	double k = dat[z][y] / dat[x][y];
	for (int i = 0; i < col; i++)
		dat[z][i] -= k * dat[x][i];
}
void Matrix::eli(int x, double div) // C[x] /= div;
{
	for (int i = 0; i < col; i++)
		dat[x][i] /= div;
}
void Matrix::Elimination(Matrix *expand)
{ // 对expand进行高斯消元。
	int f, i = 0, j = 0;
	while (i < expand->row && j < expand->col)
	{
		f = i;
		for (int k = i + 1; k < expand->row; k++)
			if (abs(expand->dat[k][j]) > abs(expand->dat[f][j]))
			{
				f = k;
			}
		for (int k = 0; k < expand->col; k++)
		{
			swap(expand->dat[i][k], expand->dat[f][k]);
		}
		if (expand->dat[i][j] == 0)
		{
			j++;
			continue;
		}
		expand->eli(i, expand->dat[i][j]); // change dat[i][j] to 1
		for (int k = i + 1; k < expand->row; ++k)
			expand->eli(i, j, k);
		i++, j++;
	}
}
void Matrix::invert(Matrix *ret, Matrix *source)
{
	int n = source->col;
	Matrix *expand = creatmatrix(n, n * 2);

	if (source->col != source->row)
		puts("Peculiar!"), exit(0);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n * 2; j++)
		{
			if (j < n)
				expand->dat[i][j] = source->dat[i][j];
			else
				expand->dat[i][j] = (j - n == i);
		}
	Elimination(expand);
	for (int i = 0; i < expand->row; i++)
		for (int j = 0; j < i; j++) // i i  C[j] -= dat[j][i] * C [i] ;
			expand->eli(i, i, j);
	for (int i = 0; i < n; i++)
		if (expand->dat[i][i] == 0)
			puts("Peculiar!"), exit(0);
	for (int i = 0; i < n; i++)
		for (int j = n; j < n * 2; j++)
			ret->dat[i][j - n] = expand->dat[i][j];
}
#define creatmatrix Matrix::creatmatrix
#define multi Matrix::multi
#define add Matrix::add
#define sub Matrix::sub
#define Trans Matrix::Trans
#define trace Matrix::trace
int n, m;
long long k;
int main(void)
{
	freopen("in", "r", stdin);
	scanf("%d%d", &n, &m);
	Matrix *A = creatmatrix(n, m), *B = creatmatrix(n, m);
	A->Input(n, m);
	B->Input(n, m);
	Matrix *C = creatmatrix(m, n);
	Trans(C, B);
	Matrix *D = creatmatrix(n, n);
	multi(D, A, C);
	D->print();
	Matrix *result = creatmatrix(n, n);
	// Matrix::quickpow(result, D, 3);
	// result->print();
	// printf("%.2lf", Matrix::Det(D));
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
		Matrix equations(p,q+1);
		printf("====== init_equations =====\n"); //���뷽����
		equations.Input(equations.row,equations.col);
		elimination(equations);
		for(int i=0;i<equations->row;++i)
			printf("x%d = %.1lf, ",i+1,Matrixat(equations,i,equations->col-1));
		return 0;
	}
	printf("�������Ľ�����");
	scanf("%d",&n);
	Matrix *A=InitMatrixatrix(n,n), *Q=InitMatrixatrix(n,n), *R=InitMatrixatrix(n,n);
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
			 AddMatrixatrix(A,Q,R);
			 printf("====== result_matrix ======\n");
			 Output(A);
	 }
	 if( *(argv[1]+1) == '-')
	 {
			 printf("====== first_matrix ======\n");
			 Input(Q,n,n);
			 printf("====== second_matrix ======\n");
			 Input(R,n,n);
			 DecMatrixatrix(A,Q,R);
			 printf("====== result_matrix ======\n");
			 Output(A);
	 }
	 if( *(argv[1]+1) == '*')
	 {
			 printf("====== first_matrix ======\n");
			 Input(Q,n,n);
			 printf("====== second_matrix ======\n");
			 Input(R,n,n);
			 MatrixulMatrixatrix(A,Q,R);
			 printf("====== result_matrix ======\n");
			 Output(A);
	 }
	 if(*(argv[1]+1) == 'd')
	 {
			printf("====== init_matrix ======\n");
			Input(A,n,n);
			printf("the determinant is %lf",detMatrixatrix(A));
	 }
	if(*(argv[1]+1) == 'i')
	{
		 printf("====== init_matrix ======\n");
		 Input(A,n,n);

		 Output(InverseMatrixatrix(A));
	}
	if(*(argv[1]+1) == 'p')
	{
		long long k;
		printf("====== init_matrix ======\n");
		Input(A,n,n);
		printf("====== mutiple_times ======\n");
		scanf("%lld",&k);
		PowerMatrixatrix(Q,A,k);
		printf("====== result_matrix =====\n");
		Output(Q);
	}
	if(*(argv[1]+1) == 't')
	{
		printf("====== init_matrix ======\n");
		Input(A,n,n);
		TransMatrixatrix(A);
		printf("====== result_matrix =====\n");
		Output(A);
	}
*/
// }