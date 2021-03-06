// MatrixVS.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "matrix.h"
#include <array>


int main()
{
	auto arrparams = std::vector<Complex>{Complex(12, 0), 
		Complex(1, 10) ,
		Complex(0, -1),
		Complex(20, 0) };

	Matrix m = Matrix(2, 2, arrparams);
	std::cout << "SOURCE (m t) =  " << m  << "\n\n";
	auto num = m.getOwnNumbers();

	for(auto i =0; i < 2; i++)
	{
		std::cout << num[i];
	}
	
	std::cout << "\ntrack = " << m.getTrack() << "\n\n";

	Matrix t = Matrix(2, 2, arrparams);
	std::cout << "+ = " << m  + t << "\n\n";
	std::cout << "- = " << m - t << "\n\n";

	m.transpose();
	
	std::cout << "* = " << m * t << "\n\n";

	std::cout << "* complex = " << m * Complex(15, 4) << "\n\n";


	std::cout << "RANDOM = " << Matrix::createRandom() << "\n\n";
	getchar();
    return 0;
}

