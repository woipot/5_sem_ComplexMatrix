#pragma once
#include <iostream>
#include <cmath>

using namespace std;

class Complex         // класс " омплексное число"
{
private:
	double re, im;      // действительна€ и мнима€ части

public:
	double real() const
	{
		return re;
	}

	double imagine() const
	{
		return im;
	}

	// конструкторы 

	Complex()
	{
		re = 0;
		im = 0;
	};

	Complex(double r)
	{
		re = r;
		im = 0;
	}

	Complex(double r, double i)
	{
		re = r;
		im = i;
	}

	Complex(const Complex &c)   // конструктор копировани€
	{
		re = c.re;
		im = c.im;
	}

	// деструктор
	~Complex()
	{
	}

	// остальные функции

	// ћодуль комплексного числа
	double abs()
	{
		return sqrt(re * re + im * im);
	}



	// оператор присваивани€
	Complex& operator = (const Complex &c)
	{
		re = c.real();
		im = c.imagine();

		return *this;
	}


	// оператор +=
	Complex& operator += (const Complex &c)
	{
		re += c.re;
		im += c.im;
		return *this;
	}

	Complex& operator -= (const Complex &c)
	{
		re -= c.re;
		im -= c.im;
		return *this;
	}

	// оператор сложени€
	Complex operator + (const Complex &c)
	{
		return Complex(re + c.re, im + c.im);
	}

	Complex operator + (double c)
	{
		return Complex(re + c, im);
	}


	// оператор вычитани€
	Complex operator - (const Complex &c)
	{
		return Complex(re - c.re, im - c.im);
	}

	Complex operator - ()
	{
		return Complex(-re , -im);
	}

	// оператор умножени€
	Complex operator * (const Complex &c)
	{
		return Complex(re * c.re - im * c.im, re * c.im + im * c.re);
	}

	Complex& operator *= (const Complex &c)
	{
		re = re * c.re - im * c.im;
		im = re * c.im + im * c.re;
		return *this;
	}


	Complex operator * (double c)
	{
		return Complex(re * c - im, re * c + im);
	}

	// оператор делени€
	Complex operator / (const Complex &c)
	{
		Complex temp;

		double r = c.re * c.re + c.im * c.im;
		temp.re = (re * c.re + im * c.im) / r;
		temp.im = (im * c.re - re * c.im) / r;

		return temp;
	}

	Complex& operator /= (const Complex &c)
	{
		Complex temp(*this);

		double r = c.re * c.re + c.im * c.im;
		re = (temp.re * c.re + temp.im * c.im) / r;
		im = (temp.im * c.re - temp.re * c.im) / r;

		return *this;
	}


	bool operator==(const Complex& obj2);

	Complex cSqrt() const;

	friend ostream & operator<< (ostream &, const Complex &);
	friend istream & operator>> (istream &, Complex &);

};

