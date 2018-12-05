/*
 * Knapsack.h
 *
 *  Created on: Sep 26, 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>
 */

#ifndef KNAPSACKDECODER_H_
#define KNAPSACKDECODER_H_

#include <vector>
#include <algorithm>
#include <numeric>
#include "Knapsack.h"

class KnapsackDecoder {
public:
	KnapsackDecoder(const Knapsack& knapsack);
	virtual ~KnapsackDecoder();

	double decode(std::vector<double>& chromosome) const;
	std::vector<bool> decodeSolution(std::vector<double>& chromosome) const;
	void correctChromosome(std::vector<double>& chromosome, const std::vector<bool>& solution) const;

	// Local Search
	std::vector<bool> twoSwap(const std::vector<bool>& solution);

protected:
	const Knapsack knapsack;

	template <typename SolutionContainer>
	double computeFitness(const SolutionContainer& solution) const;
	bool isValid(const std::vector<bool>& solution);
	void adjustSolution(std::vector<bool>& solution) const;

};

KnapsackDecoder::KnapsackDecoder(const Knapsack& _knapsack) : knapsack(_knapsack) {
}

KnapsackDecoder::~KnapsackDecoder() {
}

double KnapsackDecoder::decode(std::vector<double>& chromosome) const {
	std::vector<bool> solution  = decodeSolution(chromosome);

	return computeFitness(solution);
}

std::vector<bool> KnapsackDecoder::decodeSolution(std::vector<double>& chromosome) const {
	std::vector<bool> solution;

	for(std::vector<double>::const_iterator itKey = chromosome.begin(); itKey != chromosome.end(); ++itKey){
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

template <typename SolutionContainer>
double KnapsackDecoder::computeFitness(const SolutionContainer& solution) const {
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

void KnapsackDecoder::adjustSolution(std::vector<bool>& solution) const {
	std::vector<unsigned> items(knapsack.getN());
	fill(solution.begin(), solution.end(), false);

	std::iota(items.begin(), items.end(), 0); // initialize vector with increasing values

	std::random_shuffle(items.begin(), items.end());

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

void KnapsackDecoder::correctChromosome(std::vector<double>& chromosome, const std::vector<bool>& solution) const {
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

std::vector<bool> KnapsackDecoder::twoSwap(const std::vector<bool>& solution){
	std::vector<bool> s = solution;
	std::vector< std::pair<unsigned, unsigned> > pairs(knapsack.getN());

	for(unsigned i = 0; i < solution.size(); ++i){
		for(unsigned j = i+1; j < solution.size(); ++j){
			pairs[i].first = i;
			pairs[i].second = j;
		}
	}

	std::random_shuffle(pairs.begin(), pairs.end());

	double bestFitness = computeFitness(solution);
	for(unsigned i = 0; i < pairs.size();){
		unsigned a = pairs[i].first;
		unsigned b = pairs[i].second;
		swap(s[a], s[b]);

		if(isValid(s)){
			double fitness = computeFitness(s);
			if(fitness > bestFitness){
				bestFitness = fitness;
				std::random_shuffle(pairs.begin(), pairs.end());
				i = 0;
			}
			else{
				swap(s[b], s[a]);
				++i;
			}
		}
		else{
			swap(s[b], s[a]);
			++i;
		}

	}

	return s;
}

#endif /* KNAPSACKDECODER_H_ */
