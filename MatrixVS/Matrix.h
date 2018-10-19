#ifndef MATRIX_H
#define MATRIX_H

#include <complex>
#include <vector>
#include "Complex.h"

#define MAXITER 30


class Matrix
{

public:
	Matrix(const Matrix &obj);
	Matrix(size_t rowsCount, size_t columnsCount);
	Matrix(size_t rowsCount, size_t columnsCount,const std::vector<Complex> &arr);

	~Matrix();

	unsigned int columnsCount() const;
	unsigned int rowsCount() const;

	Complex getCase(const unsigned int &i, const unsigned int &j) const;
	void setCase(Complex num, unsigned int i, unsigned int j);

	Complex getDeterminant() const;

	Matrix extractMatrix(unsigned int i, unsigned int j) const;
	Matrix inverseMatrix() const;
	void transpose();
	Complex getTrack() const;
	
	Matrix getComplexConjugate() const;
	Matrix getErmiteConjugate() const;

	Complex *getOwnNumbers() const;
	Complex **getOwnVectors() const;

	bool isRectangle() const;
	bool isEqualSizes(const Matrix &obj) const;


	Matrix matrixPow(unsigned int degree) const;

	Matrix operator += (const Matrix &obj);
	Matrix operator -= (const Matrix &obj);
	Matrix operator *= (const Matrix &obj);
	Matrix operator /= (const Matrix &obj);
	Matrix operator = (const Matrix &obj);

	Matrix operator += (const Complex &num);
	Matrix operator -= (const Complex &num);
	Matrix operator *= (const Complex &num);
	Matrix operator /= (const Complex &num);
	Matrix operator =  (const Complex &num);

	static Matrix createRandom();

	friend std::ostream& operator<<(std::ostream& os, const Matrix& obj);

private:
	static Complex **createArr(size_t rowsCount, size_t columsCount, Complex defaultParam);
	static void setDefaultParams(size_t rowsCount, size_t columsCount, Complex defaultParam, Complex **obj);
	void resize(size_t newRowsCount, size_t newColumsCount);
	Complex* getDiagonale() const;
	void getNumbersAndVectors(Complex ***vectors, Complex **numbers) const;

private:
	size_t _columnsCount;
	size_t _rowsCount;

	Complex **_matrix;

};

Matrix operator + (const Matrix &obj1, const Matrix &obj2);
Matrix operator - (const Matrix &obj1, const Matrix &obj2);
Matrix operator * (const Matrix &obj1, const Matrix &obj2);
Matrix operator / (const Matrix &obj1, const Matrix &obj2);

Matrix operator + (const Matrix &obj, const Complex &num);
Matrix operator - (const Matrix &obj, const Complex &num);
Matrix operator * (const Matrix &obj, const Complex &num);
Matrix operator / (const Matrix &obj, const Complex &num);


#endif // MATRIX_H
