#pragma once
#include <iostream>
#include <vector>
#include "matrix.hpp"
using std::cout; using std::endl; using std::cin;
using std::string;
using ah::Matrix; using ah::readMatrix;

struct HMM {
	Matrix<double> A = Matrix <double>(1,1);
	Matrix<double> B = Matrix <double>(1,1);
	Matrix<double> pi = Matrix <double>(1,1);
	int numStates;
};

// Declarations.
void forward(HMM& hmm, std::vector<int>& obs, Matrix<double>& alpha);
// End declarations.

int main(void) {
	HMM hmm;
	hmm.A = readMatrix<double>();
	hmm.B = readMatrix<double>();
	hmm.pi = readMatrix<double>();
	int olength;
	hmm.numStates = hmm.A.cols();
	cin >> olength;
	std::vector<int> obs(olength);
	for (int i = 0; i < olength; ++i) {
		cin >> obs[i];
	}
	Matrix<double> alpha(olength,hmm.numStates);
	forward(hmm, obs, alpha);
	double res = 0;
	for(int state = 0; state < hmm.numStates; ++state) {
		res += alpha(olength-1,state);
	}
	cout << res << endl;
	return 0;
}

void forward(HMM& hmm, std::vector<int>& obs, Matrix<double>& alpha) {
	for(int state = 0; state < hmm.numStates; state++) {
		alpha(0,state) = hmm.pi(state) * hmm.B(state, obs[0]);
	}
	for(unsigned t = 1; t < obs.size(); t++) {
		for(int state = 0; state < hmm.numStates; state++) {
			alpha(t,state) = 0; //Just to be sure it's 0.
			for(int oldState = 0; oldState < hmm.numStates; oldState++) {
				alpha(t,state) += alpha(t-1,oldState) * hmm.A(oldState,state);
			}
			alpha(t,state) *= hmm.B(state, obs[t]);
		}
	}
}