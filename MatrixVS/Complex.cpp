#include "stdafx.h"
#include "Complex.h"
#include <complex>


// перегрузка оператора <<
ostream& operator<< (ostream &out, const Complex &c)
{
	out << "(" << c.re << ", " << c.im << ")";
	return out;
}

// перегрузка оператора >>
istream& operator>> (istream &in, Complex &c)
{
	in >> c.re >> c.im;
	return in;
}

bool Complex::operator==(const Complex& obj2)
{
	return real() == obj2.real() && imagine() == obj2.imagine();
}

Complex Complex::cSqrt() const
{
	auto phi = atan2(im, re); //+PI;
	auto r = sqrt(re*re + im*im);
	auto R = sqrt(r);
	auto Phi = (1 / 2) * phi;
	auto re = R * cos(Phi);
	auto im = R * sin(Phi);
	return Complex(re, im);
}
