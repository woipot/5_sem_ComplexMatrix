#include "stdafx.h"
#include "matrix.h"
#include <cstdarg>

Matrix::Matrix(const Matrix& obj)
{
	_columnsCount = obj._columnsCount;
	_rowsCount = obj._rowsCount;
	_matrix = createArr(_rowsCount, _columnsCount);

	for (unsigned int i = 0; i < _rowsCount; i++)
		for (unsigned int j = 0; j < _columnsCount; j++)
			_matrix[i][j] = obj._matrix[i][j];
}


Matrix::Matrix(size_t rowsCount, size_t columnsCount, ...)
{
	_columnsCount = columnsCount;
	_rowsCount = rowsCount;
	_matrix = createArr(rowsCount, columnsCount);


	const int paramsCount = columnsCount * rowsCount;

	va_list uk_arg;
	va_start(uk_arg, columnsCount);  /*  установка указателя uk_arg на  */


	int currentRow = 0;
	int currentColumn = 0;
	for (int counter = 0; counter < paramsCount; counter++)
	{
		if (currentColumn == columnsCount)
			currentColumn = 0;
		if (counter != 0 && counter % columnsCount == 0)
			currentRow++;


		auto param = va_arg(uk_arg, double);
		if(param == -1)
			return;
			

		_matrix[currentRow][currentColumn] = std::complex<double>(param);

		currentColumn++;
	}

	va_end(uk_arg);       

}

unsigned int Matrix::columnsCount() const
{
	return _columnsCount;
}

unsigned int Matrix::rowsCount() const
{
	return _rowsCount;
}

std::complex<double> Matrix::getCase(const unsigned& i, const unsigned& j) const
{
	if (i > rowsCount() || j > columnsCount())
		throw std::out_of_range("uncorrect index");

	return _matrix[i][j];
}

void Matrix::setCase(std::complex<double> num, unsigned i, unsigned j)
{
	this->_matrix[i][j] = num;
}

std::complex<double> Matrix::determinant() const
{
	auto d = std::complex<double>(0, 0);
	auto k = std::complex<double>(1, 0);


	if (!isRectangle())
		throw std::exception("matrix hasnt rectangle params");

	if (rowsCount() < 1)
		throw std::exception("Very small matrix");
	

	if (rowsCount() == 1) {
		d = _matrix[0][0];
		return d;
	}
	if (rowsCount() == 2) {
		d = _matrix[0][0] * _matrix[1][1] - _matrix[1][0] * _matrix[0][1];
		return d;
	}
	if (rowsCount() > 2) {
		for (unsigned int i = 0; i < rowsCount(); i++) {
			Matrix p = extractMatrix(i, 0);
			d = d + k * _matrix[i][0] * p.determinant();
			k = -k;
		}
	}
	return d;
}

Matrix Matrix::extractMatrix(unsigned i, unsigned j) const
{
	auto result = Matrix(rowsCount() - 1, rowsCount() - 1);
	bool di = false, dj = false;
	for (unsigned int ki = 0; ki < rowsCount() - 1; ki++) { // проверка индекса строки
		if (ki == i) di = true;
		dj = false;
		for (unsigned int kj = 0; kj < rowsCount() - 1; kj++) { // проверка индекса столбца
			if (kj == j) dj = true;
			result._matrix[ki][kj] = _matrix[ki + di][kj + dj];
		}
	}
	return result;
}

Matrix Matrix::inverseMatrix() const
{
	auto num = this->determinant();
	if (num.real() != 0) {
		Matrix mainMat(rowsCount(), columnsCount());
		for (unsigned int i = 0; i < rowsCount(); i++)
			for (unsigned int j = 0; j < columnsCount(); j++) {
				Matrix minor = extractMatrix(i, j);
				const auto setter = pow(-1.0, i + j) * minor.determinant() / num;
				
				mainMat.setCase(setter, i, j);
			}
		mainMat.transponire();
		return mainMat;
	}
	else throw std::exception("determianant = 0");
}

void Matrix::transponire()
{
	const Matrix a(*this);
	resize(columnsCount(), rowsCount());
	for (unsigned int i = 0; i < a.rowsCount(); i++)
		for (unsigned int j = 0; j < a.columnsCount(); j++)
			this->_matrix[j][i] = a._matrix[i][j];
}

bool Matrix::isRectangle() const
{
	return columnsCount() == rowsCount();
}

bool Matrix::isEqualSizes(const Matrix& obj) const
{
	return obj.rowsCount() == rowsCount() && obj.columnsCount() == columnsCount();
}

/////////////////////////////////////////////-----------------operatoers 
Matrix Matrix::operator+=(const Matrix& obj)
{
	if (!isEqualSizes(obj))
		throw std::exception("matrixes of different sizes");

	for (unsigned int i = 0; i < rowsCount(); i++)
		for (unsigned int j = 0; j < columnsCount(); j++) {
			this->_matrix[i][j] += obj._matrix[i][j];
		}
	return *this;
}

Matrix Matrix::operator-=(const Matrix& obj)
{
	if (!isEqualSizes(obj))
		throw std::exception("matrixes of different sizes");

	for (unsigned int i = 0; i < rowsCount(); i++)
		for (unsigned int j = 0; j < columnsCount(); j++) {
			this->_matrix[i][j] -= obj._matrix[i][j];
		}
	return *this;
}

Matrix Matrix::operator*=(const Matrix& obj)
{
	if (obj.rowsCount() != columnsCount())
		throw std::exception("matrixes of different sizes");

	Matrix copy(*this);
	for (unsigned int i = 0; i < rowsCount(); i++)
		for (unsigned int j = 0; j < columnsCount(); j++) {
			for (unsigned int q = 0; q < rowsCount(); q++)
				this->_matrix[i][j] = this->_matrix[i][j] + copy._matrix[i][q] * obj._matrix[q][j];
		}
	return *this;
}

Matrix Matrix::operator/=(const Matrix& obj)
{
	if (!isEqualSizes(obj))
		throw std::exception("matrixes of different sizes");

	try
	{
		return *this *= obj.inverseMatrix();
	}
	catch (std::exception err)
	{
		throw err;
	}
}

Matrix Matrix::operator=(const Matrix& obj)
{
	if (!isEqualSizes(obj))
		throw std::exception("matrixes of different sizes");

	for (unsigned int i = 0; i < rowsCount(); i++)
		for (unsigned int j = 0; j < columnsCount(); j++)
			this->_matrix[i][j] = obj._matrix[i][j];

	return *this;
}

Matrix Matrix::operator+=(const std::complex<double>& num)
{
	for (unsigned int i = 0; i < rowsCount(); i++)
		for (unsigned int j = 0; j < columnsCount(); j++) {
			this->_matrix[i][j] += num;
		}
	return *this;
}

Matrix Matrix::operator-=(const std::complex<double>& num)
{
	for (unsigned int i = 0; i < rowsCount(); i++)
		for (unsigned int j = 0; j < columnsCount(); j++) {
			this->_matrix[i][j] -= num;
		}
	return *this;
}

Matrix Matrix::operator*=(const std::complex<double>& num)
{
	for (unsigned int i = 0; i < rowsCount(); i++)
		for (unsigned int j = 0; j < columnsCount(); j++) {
			this->_matrix[i][j] *= num;
		}
	return *this;
}

Matrix Matrix::operator/=(const std::complex<double>& num)
{
	for (unsigned int i = 0; i < rowsCount(); i++)
		for (unsigned int j = 0; j < columnsCount(); j++) {
			this->_matrix[i][j] /= num;
		}
	return *this;
}

Matrix Matrix::operator=(const std::complex<double>& num)
{
	for (unsigned int i = 0; i < rowsCount(); i++)
		for (unsigned int j = 0; j < columnsCount(); j++) {
			this->_matrix[i][j] = num;
		}
	return *this;
}


/////------ simple operations
Matrix operator+(const Matrix& obj1, const Matrix& obj2)
{
	Matrix copy(obj1);

	try
	{
		return copy += obj2;
	}
	catch (std::exception err)
	{
		throw;
	}
}

Matrix operator-(const Matrix& obj1, const Matrix& obj2)
{
	Matrix copy(obj1);

	try
	{
		return copy -= obj2;
	}
	catch (std::exception err)
	{
		throw;
	}
}

Matrix operator*(const Matrix& obj1, const Matrix& obj2)
{
	Matrix copy(obj1);

	try
	{
		return copy *= obj2;
	}
	catch (std::exception err)
	{
		throw;
	}
}

Matrix operator/(const Matrix& obj1, const Matrix& obj2)
{
	Matrix copy(obj1);

	try
	{
		return copy /= obj2;
	}
	catch (std::exception err)
	{
		throw;
	}
}

Matrix operator+(const Matrix& obj, const std::complex<double>& num)
{
	Matrix copy(obj);

	return copy += num;
}

Matrix operator-(const Matrix& obj, const std::complex<double>& num)
{
	Matrix copy(obj);

	return copy -= num;
}

Matrix operator*(const Matrix& obj, const std::complex<double>& num)
{
	Matrix copy(obj);

	return copy *= num;
}

Matrix operator/(const Matrix& obj, const std::complex<double>& num)
{
	Matrix copy(obj);

	return copy /= num;
}

///////////////////////////////////////// ----------------------------------
void Matrix::resize(size_t newRowsCount, size_t newColumsCount)
{
	delete _matrix;

	_columnsCount = newColumsCount;
	_rowsCount = newRowsCount;

	_matrix = createArr(newRowsCount, newColumsCount);

}

std::complex<double> ** Matrix::createArr(size_t rowsCount, size_t columsCount)
{
	std::complex<double> **resultArr = new std::complex<double> *[rowsCount];
	for (int i = 0; i < rowsCount; i++)
	{
		resultArr[i] = new std::complex<double>[columsCount];
	}

	setDefaultParams(rowsCount, columsCount, 0, resultArr);
	return resultArr;
}

void Matrix::setDefaultParams(size_t rowsCount, size_t columsCount, double defaultParam, std::complex<double> **obj)
{
	for (int i = 0; i < rowsCount; i++)
	{
		for (int j = 0; j < columsCount; j++)
		{
			obj[i][j] = std::complex<double>(defaultParam, 0);
		}
	}
}


std::ostream& operator<<(std::ostream& os, const Matrix& obj)
{

	os << "[" << obj.rowsCount() << "," << obj.columnsCount() << "]:\n";
	for (unsigned int i = 0; i < obj.rowsCount(); i++) {
		for (unsigned int j = 0; j < obj.columnsCount(); j++)
		{
			os << obj._matrix[i][j];
			if (j < obj.rowsCount() - 1)
				os << ' ';
		}
		os << "\n";
	}
	return os;
}
