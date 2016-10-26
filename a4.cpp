#include <iostream>
#include <vector>
#include <float.h>
#include <cmath>
#include "matrix.hpp"
using std::cout; using std::endl; using std::cin;
using std::string;

struct HMM {
	Matrix<double> A = Matrix<double>(1,1);
	Matrix<double> B = Matrix<double>(1,1);
	Matrix<double> pi = Matrix<double>(1,1);
	int numStates;
};

struct Greeks {
	Matrix<double> alpha = Matrix<double>(1,1);
	Matrix<double> beta = Matrix<double>(1,1);
	Matrix<double> gamma = Matrix<double>(1,1);
	std::vector<Matrix<double>> digamma = std::vector<Matrix<double>>(1);
	std::vector<double> scale = std::vector<double>(1);
	double logProb = -DBL_MAX;
	int maxIters = 100;
	int iters = 0;
};

// Declarations.
void forward(HMM& hmm, Greeks& greeks, std::vector<int>& obs);
void backward(HMM& hmm, Greeks& greeks, std::vector<int>& obs);
void gammas(HMM& hmm, Greeks& greeks, std::vector<int>& obs);
void reestimate(HMM& hmm, Greeks& greeks, std::vector<int>& obs);
bool logTest(Greeks& greeks, std::vector<int>& obs);
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
	Greeks greeks;
	greeks.alpha = Matrix<double>(hmm.numStates, olength);
	greeks.beta = Matrix<double>(hmm.numStates, olength);
	greeks.gamma = Matrix<double>(hmm.numStates, olength);
	greeks.digamma = std::vector<Matrix<double>>(olength);
	for(int t = 0; t < olength; t++) {
		greeks.digamma[t] = Matrix<double>(hmm.numStates, hmm.numStates);
	}
	greeks.scale = std::vector<double>(olength);
	//Actual algorithm
	bool again = true;
	while(again) {
		forward(hmm,greeks,obs);
		backward(hmm,greeks,obs);
		gammas(hmm,greeks,obs);
		reestimate(hmm,greeks,obs);
		again = logTest(greeks,obs);
	}
	hmm.A.print();
	hmm.B.print();
	return 0;
}

void forward(HMM& hmm, Greeks& greeks, std::vector<int>& obs) {
	greeks.scale[0] = 0;
	for(int state = 0; state < hmm.numStates; state++) {
		greeks.alpha(state,0) = hmm.pi(state) * hmm.B(state,obs[0]);
		greeks.scale[0] += greeks.alpha(state,0);
	}
	for(int state = 0; state < hmm.numStates; state++) { //scale
		greeks.alpha(state,0) /= greeks.scale[0];
	}
	for(unsigned t = 1; t < obs.size(); t++) {
		greeks.scale[t] = 0;
		for(int state = 0; state < hmm.numStates; state++) {
			greeks.alpha(state,t) = 0; //Reset form previous run.
			for(int oldState = 0; oldState < hmm.numStates; oldState++) {
				greeks.alpha(state,t) += greeks.alpha(oldState,t-1) * hmm.A(oldState,state) * hmm.B(state, obs[t]);
				greeks.scale[t] += greeks.alpha(state,t);
			}
		}
		for(int state = 0; state < hmm.numStates; state++) { //scale
			greeks.alpha(state,t) /= greeks.scale[t];
		}
	}
}

void backward(HMM& hmm, Greeks& greeks, std::vector<int>& obs) {
	for(int state = 0; state < hmm.numStates; state++) {
		greeks.beta(state,obs.size()-1) = 1 / greeks.scale[obs.size()-1];
	}
	for(int t = obs.size()-2; t >= 0; t--) {
		for(int state = 0; state < hmm.numStates; state++) {
			greeks.beta(state,t) = 0;
			for(int nextState = 0; nextState < hmm.numStates; nextState++) {
				greeks.beta(state,t) += hmm.A(state,nextState) * hmm.B(nextState, obs[t+1]) * greeks.beta(nextState,t+1);
			}
			greeks.beta(state,t) /= greeks.scale[t]; //scale
		}
	}
}

void gammas(HMM& hmm, Greeks& greeks, std::vector<int>& obs) {
	double denom;
	for(unsigned t = 0; t < obs.size()-1; t++) {
		denom = 0;
		for(int state = 0; state < hmm.numStates; state++) {
			for(int nextState = 0; nextState < hmm.numStates; nextState++) {
				denom += greeks.alpha(state,t) * hmm.A(state,nextState) * hmm.B(nextState,obs[t+1]) * greeks.beta(nextState,t+1);
			}
		}
		for(int state = 0; state < hmm.numStates; state++) {
			greeks.gamma(state,t) = 0;
			for(int nextState = 0; nextState < hmm.numStates; nextState++) {
				greeks.digamma[t](state,nextState) = (greeks.alpha(state,t) * hmm.A(state,nextState) * hmm.B(nextState,obs[t+1]) * greeks.beta(nextState,t+1)) / denom;
				greeks.gamma(state,t) += greeks.digamma[t](state,nextState);
			}
		}
	}

	denom = 0;
	for(int state = 0; state < hmm.numStates; state++) {
		denom += greeks.alpha(state,obs.size()-1);
	}
	for(int state = 0; state < hmm.numStates; state++) {
		greeks.gamma(state,obs.size()-1) = greeks.alpha(state, obs.size()-1) / denom;
	}
}

void reestimate(HMM& hmm, Greeks& greeks, std::vector<int>& obs) {
	double denom, numer;
	for(int state = 0; state < hmm.numStates; state++) { //Reestimate pi
		hmm.pi(state) = greeks.gamma(state,0);
	}
	for(int state = 0; state < hmm.numStates; state++) { //Reestimate A
		for(int nextState = 0; nextState < hmm.numStates; nextState++) {
			numer = 0;
			denom = 0;
			for(unsigned t = 0; t < obs.size()-1; t++) {
				numer += greeks.digamma[t](state,nextState);
				denom += greeks.gamma(state,t);
			}
			hmm.A(state,nextState) = numer / denom;
		}
	}
	for(int state = 0; state < hmm.numStates; state++) { //Reestimate B
		for(unsigned o = 0; o < hmm.B.cols(); o++) {
			numer = 0;
			denom = 0;
			for(unsigned t = 0; t < obs.size(); t++) {
				if(obs[t] == (int)o) {
					numer += greeks.gamma(state,t);
				}
				denom += greeks.gamma(state,t);
			}
			hmm.B(state,o) = numer / denom;
		}
	}
}

bool logTest(Greeks& greeks, std::vector<int>& obs) {
	double logProb = 0;
	for(unsigned t = 0; t < obs.size()-1; t++) {
		logProb += log(1 / greeks.scale[t]);
	}
	logProb = -logProb;
	if(greeks.iters < greeks.maxIters && logProb > greeks.logProb) {
		greeks.iters++;
		greeks.logProb = logProb;
		return true;
	}
	return false;
}