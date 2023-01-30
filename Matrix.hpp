#ifndef Matrix_hpp
#define Matrix_hpp

const int MODE = 1e9 + 7;

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

public:
	int row, col;

	Matrix(int a, int b);
	Matrix();
	~Matrix();
	Matrix Trans();
	double tr();
	void Eig();
	Matrix Inv();
	static void multi(Matrix *C, const Matrix *matrix_A, const Matrix *matrix_B);
	static void sub(Matrix *C, const Matrix *matrix_A, const Matrix *matrix_B);
	static void add(Matrix *C, const Matrix *matrix_A, const Matrix *matrix_B);
	static void quickpow(Matrix *C, const Matrix *Base, long long exp);
	friend void Schmidt(Matrix &Q, Matrix &R, const Matrix matrix);
	void eli(int x, int y, int z);
	void eli(int x, double div);
	friend Matrix Elimination(Matrix Source);
	double Det();
	void print();
	void Input(int row, int col);
};

// Matrix operator*(const Matrix &matrix_A, const Matrix &matrix_B);

// Matrix operator-(const Matrix &matrix_A, const Matrix &matrix_B);

// Matrix operator+(const Matrix &matrix_A, const Matrix &matrix_B);

Matrix operator^(Matrix Base, long long exp);

void Schmidt(Matrix &Q, Matrix &R, const Matrix matrix);

Matrix Elimination(Matrix expand);

#endif