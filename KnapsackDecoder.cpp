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

KnapsackDecoder::KnapsackDecoder(const Knapsack& _knapsack) : knapsack(_knapsack) {
	// TODO Auto-generated constructor stub
}

KnapsackDecoder::~KnapsackDecoder() {
	// TODO Auto-generated destructor stub
}

double KnapsackDecoder::decode(vector<double>& chromosome) const {
	vector<bool> selection  = getSolution(chromosome);

	double value_sum = 0.0;
	for(unsigned i = 0; i < selection.size(); ++i){
		if(selection[i]){
			value_sum += knapsack.getValue(i);
		}
	}

	return value_sum;
}

vector<bool> KnapsackDecoder::getSolution(std::vector<double>& chromosome) const {
	vector<bool> selection;

	for(vector<double>::const_iterator itKey = chromosome.begin(); itKey != chromosome.end(); ++itKey){
		selection.push_back( *itKey > 0.5 );
	}

	double weight_sum = 0.0;
	for(unsigned i = 0; i < selection.size(); ++i){
		if(selection[i]){
			weight_sum += knapsack.getWeight(i);
		}
	}

	if(weight_sum > knapsack.getW()){
		adjustSelection(selection);
		correctChromosome(chromosome, selection);
	}

	return selection;
}

void KnapsackDecoder::adjustSelection(vector<bool>& selection) const {
	vector<unsigned> items(knapsack.getN());
	fill(selection.begin(), selection.end(), false);

	iota(items.begin(), items.end(), 0); // initialize vector with increasing values

	random_shuffle(items.begin(), items.end());
	double weight_sum = 0.0;
	for(unsigned i = 0; i < items.size(); ++i){
		if(weight_sum + knapsack.getWeight(i) < knapsack.getW()){
			weight_sum += knapsack.getWeight(i);
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
