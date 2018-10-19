#include "stdafx.h"
#include "matrix.h"
#include <cstdarg>
#include <ctime>
#include <algorithm>

Matrix::Matrix(const Matrix& obj)
{
	_columnsCount = obj._columnsCount;
	_rowsCount = obj._rowsCount;
	_matrix = createArr(_rowsCount, _columnsCount, Complex(0));

	for (unsigned int i = 0; i < _rowsCount; i++)
		for (unsigned int j = 0; j < _columnsCount; j++)
			_matrix[i][j] = obj._matrix[i][j];
}

Matrix::Matrix(size_t rowsCount, size_t columnsCount)
{
	_columnsCount = columnsCount;
	_rowsCount = rowsCount;
	_matrix = createArr(rowsCount, columnsCount, Complex(0));
}

Complex* Matrix::getOwnNumbers() const
{
	Complex **vectors, *numbers;
	getNumbersAndVectors(&vectors, &numbers);

	return numbers;
}

Complex** Matrix::getOwnVectors() const
{
	Complex **vectors, *numbers;
	getNumbersAndVectors(&vectors, &numbers);

	return vectors;
}

Matrix::Matrix(size_t rowsCount, size_t columnsCount,const std::vector<Complex> &arr)
{
	_columnsCount = columnsCount;
	_rowsCount = rowsCount;
	_matrix = createArr(rowsCount, columnsCount, Complex(0));


	const int paramsCount = columnsCount * rowsCount;

	int currentRow = 0;
	int currentColumn = 0;
	for (int counter = 0; counter < arr.size(); counter++)
	{
		if (currentColumn == columnsCount)
			currentColumn = 0;
		if (counter != 0 && counter % columnsCount == 0)
			currentRow++;


		auto param = arr.at(counter);

		_matrix[currentRow][currentColumn] = param;

		currentColumn++;
	}

}

Matrix::~Matrix()
{
	delete _matrix;
}



unsigned int Matrix::columnsCount() const
{
	return _columnsCount;
}

unsigned int Matrix::rowsCount() const
{
	return _rowsCount;
}

Complex Matrix::getCase(const unsigned& i, const unsigned& j) const
{
	if (i > rowsCount() || j > columnsCount())
		throw std::out_of_range("uncorrect index");

	return _matrix[i][j];
}

void Matrix::setCase(Complex num, unsigned i, unsigned j)
{
	this->_matrix[i][j] = num;
}



Complex Matrix::getDeterminant() const
{
	auto d = Complex(0, 0);
	auto k = Complex(1, 0);


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
			d = d + k * _matrix[i][0] * p.getDeterminant();
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
	auto num = this->getDeterminant();
	if (num.real() != 0) {
		Matrix mainMat(rowsCount(), columnsCount());
		for (unsigned int i = 0; i < rowsCount(); i++)
			for (unsigned int j = 0; j < columnsCount(); j++) {
				Matrix minor = extractMatrix(i, j);
				Complex setter = Complex(pow(-1.0, i + j)) * minor.getDeterminant() / num;
				
				mainMat.setCase(setter, i, j);
			}
		mainMat.transpose();
		return mainMat;
	}
	else throw std::exception("determianant = 0");
}

void Matrix::transpose()
{
	const Matrix a(*this);
	resize(columnsCount(), rowsCount());
	for (unsigned int i = 0; i < a.rowsCount(); i++)
		for (unsigned int j = 0; j < a.columnsCount(); j++)
			this->_matrix[j][i] = a._matrix[i][j];
}

Complex Matrix::getTrack() const
{
	auto result = Complex(0, 0);

	auto i = 0;
	while (i < rowsCount() && i < columnsCount())
	{
		result += _matrix[i][i];
		i++;
	}
	return result;

}


Matrix Matrix::getComplexConjugate() const
{
	Matrix copy(*this);
	auto additional = Complex(0, -1);

	for (unsigned int i = 0; i < rowsCount(); i++)
		for (unsigned int j = 0; j < columnsCount(); j++)
			copy._matrix[i][j] = copy._matrix[i][j] * additional;

	return copy;
}

Matrix Matrix::getErmiteConjugate() const
{
	auto resultMatrix = getComplexConjugate();
	resultMatrix.transpose();
	return resultMatrix;
}


bool Matrix::isRectangle() const
{
	return columnsCount() == rowsCount();
}

bool Matrix::isEqualSizes(const Matrix& obj) const
{
	return obj.rowsCount() == rowsCount() && obj.columnsCount() == columnsCount();
}


Matrix Matrix::matrixPow(unsigned degree) const
{
	if(!isRectangle())
		throw std::exception("matrix is not rectangle");
	if(degree == 0)
		throw std::exception("zero degree");

	Matrix copy(*this);
	for (auto i = 0; i < degree; i++)
	{
		copy *= copy;
	}
	
	return copy;
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

Matrix Matrix::operator+=(const Complex& num)
{
	for (unsigned int i = 0; i < rowsCount(); i++)
		for (unsigned int j = 0; j < columnsCount(); j++) {
			this->_matrix[i][j] += num;
		}
	return *this;
}

Matrix Matrix::operator-=(const Complex& num)
{
	for (unsigned int i = 0; i < rowsCount(); i++)
		for (unsigned int j = 0; j < columnsCount(); j++) {
			this->_matrix[i][j] -= num;
		}
	return *this;
}

Matrix Matrix::operator*=(const Complex& num)
{
	for (unsigned int i = 0; i < rowsCount(); i++)
		for (unsigned int j = 0; j < columnsCount(); j++) {
			this->_matrix[i][j] *= num;
		}
	return *this;
}

Matrix Matrix::operator/=(const Complex& num)
{
	for (unsigned int i = 0; i < rowsCount(); i++)
		for (unsigned int j = 0; j < columnsCount(); j++) {
			this->_matrix[i][j] /= num;
		}
	return *this;
}

Matrix Matrix::operator=(const Complex& num)
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

Matrix operator+(const Matrix& obj, const Complex& num)
{
	Matrix copy(obj);

	return copy += num;
}

Matrix operator-(const Matrix& obj, const Complex& num)
{
	Matrix copy(obj);

	return copy -= num;
}

Matrix operator*(const Matrix& obj, const Complex& num)
{
	Matrix copy(obj);

	return copy *= num;
}

Matrix operator/(const Matrix& obj, const Complex& num)
{
	Matrix copy(obj);

	return copy /= num;
}

///////////////////////////////////////// ----------------------------------
Matrix Matrix::createRandom()
{
	srand(time(nullptr));
	auto rowsCount = rand() % 50 + 1;
	auto columsCount = rand() % 50 + 1;

	auto paramsCount = rowsCount * columsCount;
	auto randParams = std::vector<Complex>();

	for (auto i = 0; i < paramsCount; i++ )
	{
		auto randomComplex = (rand(), rand());
		randParams.push_back(randomComplex);
	}


	return Matrix(rowsCount, columsCount, randParams);

}


void Matrix::resize(size_t newRowsCount, size_t newColumsCount)
{
	delete _matrix;

	_columnsCount = newColumsCount;
	_rowsCount = newRowsCount;

	_matrix = createArr(newRowsCount, newColumsCount, Complex(0));

}

Complex ** Matrix::createArr(size_t rowsCount, size_t columsCount, Complex defaultParam)
{
	Complex **resultArr = new Complex *[rowsCount];
	for (int i = 0; i < rowsCount; i++)
	{
		resultArr[i] = new Complex[columsCount];
	}

	setDefaultParams(rowsCount, columsCount, defaultParam, resultArr);
	return resultArr;
}

void Matrix::setDefaultParams(size_t rowsCount, size_t columsCount, Complex defaultParam, Complex **obj)
{
	for (int i = 0; i < rowsCount; i++)
	{
		for (int j = 0; j < columsCount; j++)
		{
			obj[i][j] = defaultParam;
		}
	}
}

Complex* Matrix::getDiagonale() const
{
	auto size = std::min(rowsCount(), columnsCount());
	auto result = new Complex[size];

	for(auto i = 0; i < size; i++)
	{
		result[i] = _matrix[i][i];
	}
	return result;
}


void Matrix::getNumbersAndVectors(Complex ***vectors, Complex **numbers) const
{
	if (!isRectangle())
		throw std::exception("Is'nt rectangle matrix");

	auto size = rowsCount();
	auto d = getDiagonale();
	auto z = createArr(size, size, Complex(1));


	auto underDCount = (size * size - size) / 2;
	Complex *e = new Complex[underDCount + 2];

	auto counter = 0;
	for (auto i = 0; i < size; i++)
	{
		for (auto j = 0; j < counter; j++)
		{
			e[i + 2] = _matrix[i][j];
		}
		counter++;
	}


	int m, l, iter, i, k;
	Complex s, r, p, g, f, dd, c, b;
	for (i = 2; i <= size; i++) e[i - 1] = e[i];
	e[size] = 0.;
	for (l = 1; l <= size; l++) {
		iter = 0;
		do {
			for (m = l; m <= size - 1; m++) {
				dd = d[m].abs() + d[m + 1].abs();
				if ((dd + e[m].abs() ) == dd)
				{
					break;
				}
			}
			if (m != l) {
				if (++iter >= MAXITER) throw std::exception("Too many iterations in tqli");
				g = (d[l + 1] - d[l]) / (e[l] * 2);


				r = (Complex(1) + g * g).cSqrt();
				if (g.real() >= 0.) g += r.abs();
				else g -= r.abs();
				g = d[m] - d[l] + e[l] / g;
				s = c = 1.; p = 0.;
				for (i = m - 1; i >= l; i--) {
					f = s * e[i]; b = c * e[i];
					e[i + 1] = r = (f*f + g * g).cSqrt();
					if (r == 0.) { d[i + 1] -= p; e[m] = 0.; break; }
					s = f / r; c = g / r; g = d[i + 1] - p; r = (d[i] - g)*s + c*b*2;
					d[i + 1] = g + (p = s * r); g = c * r - b;
					for (k = 1; k <= size; k++) {
						f = z[k][i + 1]; z[k][i + 1] = s * z[k][i] + c * f; z[k][i] = c * z[k][i] - s * f;
					}
				}
				if (r == 0. && i >= l) continue;
				d[l] -= p; e[l] = g; e[m] = 0.;
			}
		} while (m != l);
	}

	*numbers = d;
	*vectors = z;
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
