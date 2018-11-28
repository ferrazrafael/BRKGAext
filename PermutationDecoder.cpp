/*
 * PermutationDecoder.cpp
 *
 *  Created on: Nov 20, 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>
 */

#include <algorithm>
#include "PermutationDecoder.h"

using namespace std;

PermutationDecoder::PermutationDecoder(const Permutation _permutation) : permutation(_permutation) {
}

PermutationDecoder::~PermutationDecoder() {
}

double PermutationDecoder::decode(std::vector<double>& chromosome) const {
	vector<ItemType> solution  = decodeSolution(chromosome);

	return computeFitness(solution);
}

double PermutationDecoder::computeFitness(const std::vector<ItemType>& solution) const {
	// TODO API user implements fitness computation based on problem data
}

vector<ItemType> PermutationDecoder::decodeSolution(vector<double>& chromosome) const {
	vector<ItemType> solution = permutation.getItems();
	vector< pair<double, ItemType> > pairs( permutation.getN() );

	for(unsigned i = 0; i < permutation.getN(); ++i){
		pairs[i].first = chromosome[i];
		pairs[i].second = solution[i];
	}

	// sort by key
	sort(pairs.begin(), pairs.end(),
			[](pair<double, ItemType> a, pair<double, ItemType> b) { return a.first < b.first; } );

	// copy items on order to solution
	transform(pairs.begin(), pairs.end(), solution.begin(),
			[](pair<double, ItemType> p) { return p.second; } );

	/*
	 * TODO Implement if the problem has restrictions
	 *
	 if(Some restriction is violated){
		adjustSolution(solution);
		correctChromosome(chromosome, solution);
	}
	 */

	return solution;
}

void PermutationDecoder::correctChromosome(vector<double>& chromosome, const vector<ItemType>& solution) const {
	vector<ItemType> items = permutation.getItems();
	vector< pair<double, ItemType> > pairs( permutation.getN() );

	for(unsigned i = 0; i < permutation.getN(); ++i){
		pairs[i].first = chromosome[i];
		pairs[i].second = items[i];
	}

	// sort by key
	sort(pairs.begin(), pairs.end(),
			[](pair<double, ItemType> a, pair<double, ItemType> b) { return a.first < b.first; } );

	for(unsigned i = 0; i < permutation.getN(); ++i){
		if(pairs[i].second != solution[i]){
			for(unsigned j = i+1; j < permutation.getN(); ++j){
				if(solution[i] == pairs[j].first){
					swap(chromosome[i], chromosome[j]);
				}
			}
		}
	}
}

vector<ItemType> PermutationDecoder::twoSwap(const vector<ItemType>& solution) {
	vector<ItemType> s = solution;

	double bestFitness = computeFitness(solution);
	for(unsigned i = 0; i < solution.size(); ++i){
		for(unsigned j = i+1; j < solution.size(); ++j){
			swap(s[i], s[j]);

			if(isValid(s)){
				double fitness = computeFitness(s);
				if(fitness > bestFitness){
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

bool PermutationDecoder::isValid(const std::vector<ItemType>& solution) {
	// TODO Implement problem requirements validation
}

vector<ItemType> PermutationDecoder::kSwap(const vector<ItemType>& solution) {
}

vector<ItemType> PermutationDecoder::bestInsertion(const vector<ItemType>& solution) {
}

vector<ItemType> PermutationDecoder::twoOpt(const vector<ItemType>& solution) {
}

void PermutationDecoder::adjustSolution(vector<ItemType>& solution) const {
	// TODO Implement to adjust invalid solution
}
