#pragma once
#include <iostream>
#include <vector>
#include "matrix.hpp"
using std::cout; using std::endl; using std::cin;
using std::string;
//using ah::Matrix; using ah::readMatrix;

struct HMM {
	Matrix<double> A = Matrix <double>(1,1);
	Matrix<double> B = Matrix <double>(1,1);
	Matrix<double> pi = Matrix <double>(1,1);
	int numStates;
};

// Declarations.
void viterbi(HMM& hmm, std::vector<int>& obs, Matrix<double>& delta, Matrix<int>& idx);
void printSolution(Matrix<double>& delta, Matrix<int>& idx, int olength, int numStates);
// End declarations.

int main(void) {
	HMM hmm;
	hmm.A = readMatrix<double>();
	hmm.B = readMatrix<double>();
	hmm.pi = readMatrix<double>();
	hmm.numStates = hmm.A.cols();
	int olength;
	cin >> olength;
	std::vector<int> obs(olength);
	for (int i = 0; i < olength; ++i) {
		cin >> obs[i];
	}
	Matrix<double> delta(olength,hmm.numStates);
	Matrix<int> idx(olength-1,hmm.numStates);
	viterbi(hmm, obs, delta, idx);
	//delta.print();
	//idx.print();
	printSolution(delta, idx, olength, hmm.numStates);
	return 0;
}

void viterbi(HMM& hmm, std::vector<int>& obs, Matrix<double>& delta, Matrix<int>& idx) {
	for(int state = 0; state < hmm.numStates; state++) {
		delta(0,state) = hmm.pi(state) * hmm.B(state, obs[0]);
	}
	double prob;
	for(unsigned t = 1; t < obs.size(); t++) {
		for(int state = 0; state < hmm.numStates; state++) {
			for(int oldState = 0; oldState < hmm.numStates; oldState++) {
				prob = hmm.A(oldState,state) * delta(t-1,oldState);
				if(prob > delta(t,state)) {
					delta(t,state) = prob;
					idx(t-1,state) = oldState;
				}
			}
			delta(t,state) *= hmm.B(state, obs[t]);
		}
	}
}

void printSolution(Matrix<double>& delta, Matrix<int>& idx, int olength, int numStates) {
	std::vector<int> solution(olength);
	double prev = -1;
	for(int state = 0; state < numStates; ++state) {
		if(delta(olength-1,state) > prev) {
			prev = delta(olength-1,state);
			solution[olength-1] = state;
		}
	}
	for(int t = olength-1; t > 0; --t) {
		int state = solution[t];
		solution[t-1] = idx(t-1, state);
	}

	/*
	for(int t = olength-2; t >= 0; t--) {
		solution[t] = idx(t,solution[t+1]);
	}*/
	for(unsigned i = 0; i < solution.size(); i++) {
		cout << solution[i] << " ";
	}
	cout << endl;
}