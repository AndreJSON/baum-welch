#pragma once
#include <iostream>
#include "matrix.hpp"
using std::cout; using std::endl; using std::cin;
using std::string;
using ah::Matrix; using ah::readMatrix;

int main(void) {
	Matrix<double> A = readMatrix<double>();
	Matrix<double> B = readMatrix<double>();
	Matrix<double> pi = readMatrix<double>();
	Matrix<double> res = pi * A * B;
	res.print();
	return 0;
}