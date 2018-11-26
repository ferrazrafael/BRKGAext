/*
 * Knapsack.cpp
 *
 *  Created on: Sep 26, 2018
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
}

KnapsackDecoder::~KnapsackDecoder() {
}

double KnapsackDecoder::decode(vector<double>& chromosome) const {
	vector<bool> solution  = decodeSolution(chromosome);

	return computeFitness(solution);
}

vector<bool> KnapsackDecoder::decodeSolution(std::vector<double>& chromosome) const {
	vector<bool> solution;

	for(vector<double>::const_iterator itKey = chromosome.begin(); itKey != chromosome.end(); ++itKey){
		solution.push_back( *itKey > 0.5 );
	}

	double weight_sum = 0.0;
	for(unsigned i = 0; i < solution.size(); ++i){
		if(solution[i]){
			weight_sum += knapsack.getWeight(i);
		}
	}

	if(weight_sum > knapsack.getW()){
		adjustSolution(solution);
		correctChromosome(chromosome, solution);
	}

	return solution;
}

double KnapsackDecoder::computeFitness(const std::vector<bool>& solution) const {
	double value_sum = 0.0;
	for(unsigned i = 0; i < solution.size(); ++i){
		if(solution[i]){
			value_sum += knapsack.getValue(i);
		}
	}

	return value_sum;
}

bool KnapsackDecoder::isValid(const std::vector<bool>& solution) {
	double weight = 0.0;
	for(unsigned i = 0; i < solution.size(); ++i){
		weight += knapsack.getWeight(i) * solution[i];
	}

	return weight <= knapsack.getW();
}

void KnapsackDecoder::adjustSolution(vector<bool>& solution) const {
	vector<unsigned> items(knapsack.getN());
	fill(solution.begin(), solution.end(), false);

	iota(items.begin(), items.end(), 0); // initialize vector with increasing values

	random_shuffle(items.begin(), items.end());
	double weight_sum = 0.0;
	for(unsigned i = 0; i < items.size(); ++i){
		if(weight_sum + knapsack.getWeight(i) < knapsack.getW()){
			weight_sum += knapsack.getWeight(i);
			solution[ items[i] ] = true;
		}
		else{
			break;
		}
	}
}

void KnapsackDecoder::correctChromosome(vector<double>& chromosome, const vector<bool>& solution) const {
	for(unsigned i = 0; i < chromosome.size(); ++i){
		if(solution[i]){
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

vector<bool> KnapsackDecoder::twoSwap(const vector<bool>& solution){
	vector<bool> s = solution;
	double value = 0.0;
	double weight = 0.0;

	double bestFitness = computeFitness(solution);
	for(unsigned i = 0; i < solution.size(); ++i){
		for(unsigned j = i+1; j < solution.size(); ++j){
			swap(s[i], s[j]);

			if(isValid(s)){ // checks if solution meets problem requirements
				double fitness = computeFitness(s);
				if(weight < knapsack.getW() && fitness > bestFitness){
					bestFitness = fitness;
				}
				else{
					swap(s[j], s[i]);
				}
			}
			else{
				swap(s[j], s[i]);
			}
		}
	}


	return s;
}
