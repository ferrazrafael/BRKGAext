/*
 * Knapsack.cpp
 *
 *  Created on: Oct 4, 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>
 */

#include "Knapsack.h"

using namespace std;

Knapsack::Knapsack(unsigned _n, const vector<double>& _values, const vector<double>& _weights, double _W) :
			n(_n), values(_values), weights(_weights), W(_W) {

}

Knapsack::~Knapsack() {
	// TODO Auto-generated destructor stub
}

double Knapsack::getValue(unsigned i) const {
	return values[i];
}

unsigned Knapsack::getN() const {
	return n;
}

double Knapsack::getW() const {
	return W;
}

const std::vector<double>& Knapsack::getValues() const {
	return values;
}

const std::vector<double>& Knapsack::getWeights() const {
	return weights;
}

double Knapsack::getWeight(unsigned i) const {
	return weights[i];
}
