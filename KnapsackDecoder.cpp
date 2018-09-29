/*
 * Knapsack.cpp
 *
 *  Created on: 26 de set de 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>
 */

#include <list>
#include <set>
#include <algorithm>
#include <numeric>
#include <random>

#include "KnapsackDecoder.h"

using namespace std;

KnapsackDecoder::KnapsackDecoder(unsigned _n, std::vector<double> _values, std::vector<double> _weights, double _W)
: n(_n), values(_values), weights(_weights), W(_W) {
	// TODO Auto-generated constructor stub
}

KnapsackDecoder::~KnapsackDecoder() {
	// TODO Auto-generated destructor stub
}

double KnapsackDecoder::decode(vector<double>& chromosome) const {
	vector<bool> selection;

	for(vector<double>::const_iterator itKey = chromosome.begin(); itKey != chromosome.end(); ++itKey){
		selection.push_back( *itKey > 0.5 );
	}

	double weight_sum = 0.0;
	double value_sum = 0.0;
	for(unsigned i = 0; i < selection.size(); ++i){
		if(selection[i]){
			weight_sum += weights[i];
		}
	}

	if(weight_sum > W){
		adjustSelection(selection);
		correctChromosome(chromosome, selection);
	}

	for(unsigned i = 0; i < selection.size(); ++i){
		if(selection[i]){
			value_sum += values[i];
		}
	}

	return value_sum;
}

unsigned KnapsackDecoder::getN() const {
	return n;
}

const std::vector<double>& KnapsackDecoder::getValues() const {
	return values;
}

double KnapsackDecoder::getW() const {
	return W;
}

const std::vector<double>& KnapsackDecoder::getWeights() const {
	return weights;
}

void KnapsackDecoder::adjustSelection(vector<bool>& selection) const {
	vector<unsigned> items(n);
	fill(selection.begin(), selection.end(), false);

	iota(items.begin(), items.end(), 0); // initialize vector with increasing values

	shuffle(items.begin(), items.end(), default_random_engine(0));
	double weight_sum = 0.0;
	for(unsigned i = 0; i < items.size(); ++i){
		if(weight_sum + weights[i] < W){
			weight_sum += weights[i];
			selection[ items[i] ] = true;
		}
		else{
			break;
		}
	}
}

void KnapsackDecoder::correctChromosome(vector<double>& chromosome, vector<bool>& selection) const {
	for(unsigned i = 0; i < chromosome.size(); ++i){
		if(selection[i]){
			if(chromosome[i] <= 0.5){
				chromosome[i] += 0.5;
			}
		}
		else{
			if(chromosome[i] > 0.5){
				chromosome[i] -= 0.5;
			}
		}
	}
}
