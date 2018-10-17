#ifndef MATRIX_H
#define MATRIX_H


#include <iostream>
#include <complex>
#include <vector>

class Matrix
{

public:
	Matrix(const Matrix &obj);
	Matrix(size_t rowsCount, size_t columnsCount, ...);

	//friend std::ostream& operator <<(std::ostream &os, const Matrix &obj);

	unsigned int columnsCount() const;
	unsigned int rowsCount() const;

	std::complex<double> getCase(const unsigned int &i, const unsigned int &j) const;
	void setCase(std::complex<double> num, unsigned int i, unsigned int j);

	std::complex<double> determinant() const;

	Matrix extractMatrix(unsigned int i, unsigned int j) const;
	Matrix inverseMatrix() const;
	void transponire();

	bool isRectangle() const;
	bool isEqualSizes(const Matrix &obj) const;


	Matrix operator += (const Matrix &obj);
	Matrix operator -= (const Matrix &obj);
	Matrix operator *= (const Matrix &obj);
	Matrix operator /= (const Matrix &obj);
	Matrix operator = (const Matrix &obj);

	Matrix operator += (const std::complex<double> &num);
	Matrix operator -= (const std::complex<double> &num);
	Matrix operator *= (const std::complex<double> &num);
	Matrix operator /= (const std::complex<double> &num);
	Matrix operator =  (const std::complex<double> &num);

	friend std::ostream& operator<<(std::ostream& os, const Matrix& obj);

private:
	static std::complex<double> **createArr(size_t rowsCount, size_t columsCount);
	static void setDefaultParams(size_t rowsCount, size_t columsCount, double defaultParam, std::complex<double> **obj);
	void resize(size_t newRowsCount, size_t newColumsCount);

private:
	size_t _columnsCount;
	size_t _rowsCount;

	std::complex<double> **_matrix;

};

Matrix operator + (const Matrix &obj1, const Matrix &obj2);
Matrix operator - (const Matrix &obj1, const Matrix &obj2);
Matrix operator * (const Matrix &obj1, const Matrix &obj2);
Matrix operator / (const Matrix &obj1, const Matrix &obj2);

Matrix operator + (const Matrix &obj, const std::complex<double> &num);
Matrix operator - (const Matrix &obj, const std::complex<double> &num);
Matrix operator * (const Matrix &obj, const std::complex<double> &num);
Matrix operator / (const Matrix &obj, const std::complex<double> &num);


#endif // MATRIX_H
