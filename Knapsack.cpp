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

vector<bool> twoSwap(const vector<bool>& selection, const Knapsack& knapsack) {
	vector<bool> s = selection;
	double value = 0.0;
	double weight = 0.0;
	const vector<double>& values = knapsack.getValues();
	const vector<double>& weights = knapsack.getWeights();

	for(unsigned i = 0; i < selection.size(); ++i){
		value += values[i] * selection[i];
		weight += weights[i] * selection[i];
	}

	double bestValue = value;
	for(unsigned i = 0; i < selection.size(); ++i){
		for(unsigned j = 1; j < selection.size(); ++j){
			swap(s[i], s[j]);
			value -= values[i] * selection[i];
			value -= values[j] * selection[j];
			weight -= weights[i] * selection[i];
			weight -= weights[j] * selection[j];
			value += values[i] * s[i];
			value += values[j] * s[j];
			weight += weights[i] * s[i];
			weight += weights[j] * s[j];
			if(weight < knapsack.getW() && value > bestValue){
				bestValue = value;
			}
			else{
				value -= values[i] * s[i];
				value -= values[j] * s[j];
				weight -= weights[i] * s[i];
				weight -= weights[j] * s[j];
				swap(s[j], s[i]);
				value += values[i] * s[i];
				value += values[j] * s[j];
				weight += weights[i] * s[i];
				weight += weights[j] * s[j];
			}
		}
	}

	return s;
}
