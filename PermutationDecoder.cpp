/*
 * PermutationDecoder.cpp
 *
 *  Created on: Nov 20, 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>
 */

#include <algorithm>
#include <list>
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
	vector< pair<unsigned, unsigned> > pairs(permutation.getN());

	for(unsigned i = 0; i < permutation.getN(); ++i){
		for(unsigned j = i+1; j < permutation.getN(); ++j){
			pairs[i].first = i;
			pairs[i].second = j;
		}
	}

	random_shuffle(pairs.begin(), pairs.end());

	double bestFitness = computeFitness(solution);
	for(unsigned i = 0; i < pairs.size();){
		unsigned a = pairs[i].first;
		unsigned b = pairs[i].second;
		swap(s[a], s[b]);

		if(isValid(s)){
			double fitness = computeFitness(s);
			if(fitness > bestFitness){
				bestFitness = fitness;
				random_shuffle(pairs.begin(), pairs.end());
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

bool PermutationDecoder::isValid(const std::vector<ItemType>& solution) {
	// TODO Implement problem requirements validation
}

vector<ItemType> PermutationDecoder::kSwap(const vector<ItemType>& solution, unsigned k) {
	vector<ItemType> s = solution;
	vector< pair<unsigned, unsigned> > pairs(permutation.getN());

	for(unsigned i = 0; i < permutation.getN(); ++i){
		for(unsigned j = i+1; j < permutation.getN(); ++j){
			pairs[i].first = i;
			pairs[i].second = j;
		}
	}

	random_shuffle(pairs.begin(), pairs.end());

	double bestFitness = computeFitness(solution);
	for(unsigned i = 0; i < pairs.size();){
		for(unsigned j = 0, c = i; j < k; ++j, ++c){
			unsigned x = c % pairs.size();
			unsigned a = pairs[x].first;
			unsigned b = pairs[x].second;
			swap(s[a], s[b]);
		}

		if(isValid(s)){
			double fitness = computeFitness(s);
			if(fitness > bestFitness){
				bestFitness = fitness;
				random_shuffle(pairs.begin(), pairs.end());
				i = 0;
			}
			else{
				for(unsigned j = k+1, c = i+j-1; j > 0; --j, --c){
					unsigned x = c % pairs.size();
					unsigned a = pairs[x].first;
					unsigned b = pairs[x].second;
					swap(s[b], s[a]);
				}
				++i;
			}
		}
		else{
			for(unsigned j = k+1, c = i+j-1; j > 0; --j, --c){
				unsigned x = c % pairs.size();
				unsigned a = pairs[x].first;
				unsigned b = pairs[x].second;
				swap(s[b], s[a]);
			}
			++i;
		}
	}

	return s;
}

vector<ItemType> PermutationDecoder::bestInsertion(const vector<ItemType>& solution) {
	list<ItemType> solutionList(solution.begin(), solution.end()); // copies solution elements to s list

	vector<ItemType> items = solution;
	random_shuffle(items.begin(), items.end());

	vector<ItemType> s = solution;
	double bestFitness = computeFitness(solution);
	for(unsigned i = 0; i < items.size(); ++i){
		list<ItemType>::iterator pos = find(solutionList.begin(), solutionList.end(), items[i]);
		solutionList.erase(pos);

		for(list<ItemType>::iterator it = solutionList.begin(); it != solutionList.end(); ++it){
			solutionList.insert(it, items[i]);
			vector<ItemType> solutionVector(solutionList.begin(), solutionList.end());
			double fitness = computeFitness(solutionVector);
			if(fitness > bestFitness){
				bestFitness = fitness;
				s = solutionVector;
			}
			solutionList.remove(items[i]);
		}

		solutionList.insert(--pos, items[i]);
	}

	return s;
}

vector<ItemType> PermutationDecoder::twoOpt(const vector<ItemType>& solution) {
	vector<ItemType> s = solution;
	vector< pair< vector<ItemType>::iterator, vector<ItemType>::iterator> > pairs;
	pairs.reserve(s.size());

	for(vector<ItemType>::iterator itA = s.begin(); itA != s.end(); ++itA){
		for(vector<ItemType>::iterator itB = next(itA); itB != s.end(); ++itB){
			pairs.push_back( make_pair(itA, itB) );
		}
	}

	random_shuffle(pairs.begin(), pairs.end());

	double bestFitness = computeFitness(solution);
	for(unsigned i = 0; i < pairs.size(); ++i){
		auto itA = pairs[i].first;
		auto itB = pairs[i].second;
		reverse(itA, itB);

		if(isValid(s)){
			double fitness = computeFitness(s);
			if(fitness > bestFitness){
				bestFitness = fitness;
			}
			else{
				reverse(itB, itA);
			}
		}
		else{
			reverse(itB, itA);
		}

	}

	return s;
}

void PermutationDecoder::adjustSolution(vector<ItemType>& solution) const {
	// TODO Implement to adjust invalid solution
}
