#pragma once
#include <iostream>
#include <cmath>

using namespace std;

class Complex         // ����� "����������� �����"
{
private:
	double re, im;      // �������������� � ������ �����

public:
	double real() const
	{
		return re;
	}

	double imagine() const
	{
		return im;
	}

	// ������������ 

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

	Complex(const Complex &c)   // ����������� �����������
	{
		re = c.re;
		im = c.im;
	}

	// ����������
	~Complex()
	{
	}

	// ��������� �������

	// ������ ������������ �����
	double abs()
	{
		return sqrt(re * re + im * im);
	}



	// �������� ������������
	Complex& operator = (const Complex &c)
	{
		re = c.real();
		im = c.imagine();

		return *this;
	}


	// �������� +=
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

	// �������� ��������
	Complex operator + (const Complex &c)
	{
		return Complex(re + c.re, im + c.im);
	}

	Complex operator + (double c)
	{
		return Complex(re + c, im);
	}


	// �������� ���������
	Complex operator - (const Complex &c)
	{
		return Complex(re - c.re, im - c.im);
	}

	Complex operator - ()
	{
		return Complex(-re , -im);
	}

	// �������� ���������
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

	// �������� �������
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

