#ifndef Matrix_hpp
#define Matrix_hpp

class Vec
{
protected: // nonstatic member
	double *dat;

public: // static member
	int row;

	Vec(int a);
	Vec();
	~Vec();
	void input();
	void print();
	double &operator[](const int x); // operator [] 必须是成员函数？
};
Vec operator+(Vec &A, Vec &B); // 全局函数
Vec operator-(Vec &A, Vec B);
double inprod(Vec &A, Vec &B);
Vec proj(Vec &vecA, Vec &vecB);
class Matrix
{
protected:
	double **dat;
	int row, col; // 非静态，每个对象专属  如果是静态，整个类共享。

public:
	static Matrix *creatmatrix(int row, int col);
	static void recycle(Matrix *matrix);
	static void Trans(Matrix *result, const Matrix *matrix);
	double tr();
	void Eig();
	static void invert(Matrix *ret, Matrix *source);
	static void multi(Matrix *C, const Matrix *matrix_A, const Matrix *matrix_B);
	static void sub(Matrix *C, const Matrix *matrix_A, const Matrix *matrix_B);
	static void add(Matrix *C, const Matrix *matrix_A, const Matrix *matrix_B);
	static void quickpow(Matrix *C, Matrix *Base, long long exp);
	static void Schmidt(Matrix *Q, Matrix *R, const Matrix *matrix);
	void eli(int x, int y, int z);
	void eli(int x, double div);
	static void Elimination(Matrix *expand); // 静态成员函数。在整个类中共享。对类内的函数可以直接使用；类外使用需要 类名::方法名
	double Det();
	void print(); // 非静态成员函数，使用时必须与一个对象匹配。如 A.print()
	void Input(int row, int col);
} M;

#endif