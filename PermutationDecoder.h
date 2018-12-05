/*
 * PermutationDecoder.h
 *
 *  Created on: Nov 20, 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>
 */

#ifndef PERMUTATIONDECODER_H_
#define PERMUTATIONDECODER_H_

#include <vector>
#include <algorithm>
#include <list>
#include "Permutation.h"

class PermutationDecoder {
public:
	PermutationDecoder(const Permutation permutation);
	virtual ~PermutationDecoder();

	double decode(std::vector<double>& chromosome) const;
	std::vector<ItemType> decodeSolution(std::vector<double>& chromosome) const;
	void correctChromosome(std::vector<double>& chromosome, const std::vector<ItemType>& solution) const;

	// Local Search
	std::vector<ItemType> twoSwap(const std::vector<ItemType>& solution);
	std::vector<ItemType> kSwap(const std::vector<ItemType>& solution, unsigned k);
	std::vector<ItemType> bestInsertion(const std::vector<ItemType>& solution);
	std::vector<ItemType> twoOpt(const std::vector<ItemType>& solution);

protected:
	const Permutation permutation;

	template <typename SolutionContainer>
	double computeFitness(const SolutionContainer& solution) const;
	bool isValid(const std::vector<ItemType>& solution);
	void adjustSolution(std::vector<ItemType>& solution) const;
};

PermutationDecoder::PermutationDecoder(const Permutation _permutation) : permutation(_permutation) {
}

PermutationDecoder::~PermutationDecoder() {
}

double PermutationDecoder::decode(std::vector<double>& chromosome) const {
	std::vector<ItemType> solution  = decodeSolution(chromosome);

	return computeFitness(solution);
}

template <typename SolutionContainer>
double PermutationDecoder::computeFitness(const SolutionContainer& solution) const {
	// TODO API user implements fitness computation based on problem data
}

std::vector<ItemType> PermutationDecoder::decodeSolution(std::vector<double>& chromosome) const {
	std::vector<ItemType> solution = permutation.getItems();
	std::vector< std::pair<double, ItemType> > pairs( permutation.getN() );

	for(unsigned i = 0; i < permutation.getN(); ++i){
		pairs[i].first = chromosome[i];
		pairs[i].second = solution[i];
	}

	// sort by key
	std::sort(pairs.begin(), pairs.end(),
			[](std::pair<double, ItemType> a, std::pair<double, ItemType> b) { return a.first < b.first; } );

	// copy items on order to solution
	transform(pairs.begin(), pairs.end(), solution.begin(),
			[](std::pair<double, ItemType> p) { return p.second; } );

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

void PermutationDecoder::correctChromosome(std::vector<double>& chromosome, const std::vector<ItemType>& solution) const {
	std::vector<ItemType> items = permutation.getItems();
	std::vector< std::pair<double, ItemType> > pairs( permutation.getN() );

	for(unsigned i = 0; i < permutation.getN(); ++i){
		pairs[i].first = chromosome[i];
		pairs[i].second = items[i];
	}

	// sort by key
	sort(pairs.begin(), pairs.end(),
			[](std::pair<double, ItemType> a, std::pair<double, ItemType> b) { return a.first < b.first; } );

	for(unsigned i = 0; i < permutation.getN(); ++i){
		if(pairs[i].second != solution[i]){
			for(unsigned j = i+1; j < permutation.getN(); ++j){
				if(solution[i] == pairs[j].first){
					std::swap(chromosome[i], chromosome[j]);
				}
			}
		}
	}
}

std::vector<ItemType> PermutationDecoder::twoSwap(const std::vector<ItemType>& solution) {
	std::vector<ItemType> s = solution;
	std::vector< std::pair<unsigned, unsigned> > pairs(permutation.getN());

	for(unsigned i = 0; i < permutation.getN(); ++i){
		for(unsigned j = i+1; j < permutation.getN(); ++j){
			pairs[i].first = i;
			pairs[i].second = j;
		}
	}

	std::random_shuffle(pairs.begin(), pairs.end());

	double bestFitness = computeFitness(solution);
	for(unsigned i = 0; i < pairs.size();){
		unsigned a = pairs[i].first;
		unsigned b = pairs[i].second;
		std::swap(s[a], s[b]);

		if(isValid(s)){
			double fitness = computeFitness(s);
			if(fitness > bestFitness){
				bestFitness = fitness;
				std::random_shuffle(pairs.begin(), pairs.end());
				i = 0;
			}
			else{
				std::swap(s[b], s[a]);
				++i;
			}
		}
		else{
			std::swap(s[b], s[a]);
			++i;
		}

	}

	return s;
}

bool PermutationDecoder::isValid(const std::vector<ItemType>& solution) {
	// TODO Implement problem requirements validation
}

std::vector<ItemType> PermutationDecoder::kSwap(const std::vector<ItemType>& solution, unsigned k) {
	std::vector<ItemType> s = solution;
	std::vector< std::pair<unsigned, unsigned> > pairs(permutation.getN());

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
			std::swap(s[a], s[b]);
		}

		if(isValid(s)){
			double fitness = computeFitness(s);
			if(fitness > bestFitness){
				bestFitness = fitness;
				std::random_shuffle(pairs.begin(), pairs.end());
				i = 0;
			}
			else{
				for(unsigned j = k+1, c = i+j-1; j > 0; --j, --c){
					unsigned x = c % pairs.size();
					unsigned a = pairs[x].first;
					unsigned b = pairs[x].second;
					std::swap(s[b], s[a]);
				}
				++i;
			}
		}
		else{
			for(unsigned j = k+1, c = i+j-1; j > 0; --j, --c){
				unsigned x = c % pairs.size();
				unsigned a = pairs[x].first;
				unsigned b = pairs[x].second;
				std::swap(s[b], s[a]);
			}
			++i;
		}
	}

	return s;
}

std::vector<ItemType> PermutationDecoder::bestInsertion(const std::vector<ItemType>& solution) {
	std::list<ItemType> sList(solution.begin(), solution.end()); // copies solution elements to s list

	std::vector<ItemType> items = solution;
	std::random_shuffle(items.begin(), items.end());

	std::vector<ItemType> s = solution;
	double bestFitness = computeFitness(solution);
	for(unsigned i = 0; i < items.size(); ++i){
		std::list<ItemType>::iterator pos = find(sList.begin(), sList.end(), items[i]);
		sList.erase(pos);

		for(std::list<ItemType>::iterator it = sList.begin(); it != sList.end(); ++it){
			sList.insert(it, items[i]);
			double fitness = computeFitness(sList);
			if(fitness > bestFitness){
				bestFitness = fitness;
				std::copy(sList.begin(), sList.begin(), s.begin());
			}
			sList.remove(items[i]);
		}

		sList.insert(--pos, items[i]);
	}

	return s;
}

std::vector<ItemType> PermutationDecoder::twoOpt(const std::vector<ItemType>& solution) {
	std::vector<ItemType> s = solution;
	std::vector< std::pair< std::vector<ItemType>::iterator, std::vector<ItemType>::iterator> > pairs;
	pairs.reserve(s.size());

	for(std::vector<ItemType>::iterator itA = s.begin(); itA != s.end(); ++itA){
		for(std::vector<ItemType>::iterator itB = next(itA); itB != s.end(); ++itB){
			pairs.push_back( make_pair(itA, itB) );
		}
	}

	std::random_shuffle(pairs.begin(), pairs.end());

	double bestFitness = computeFitness(solution);
	for(unsigned i = 0; i < pairs.size(); ++i){
		auto itA = pairs[i].first;
		auto itB = pairs[i].second;
		std::reverse(itA, itB);

		if(isValid(s)){
			double fitness = computeFitness(s);
			if(fitness > bestFitness){
				bestFitness = fitness;
			}
			else{
				std::reverse(itB, itA);
			}
		}
		else{
			std::reverse(itB, itA);
		}

	}

	return s;
}

void PermutationDecoder::adjustSolution(std::vector<ItemType>& solution) const {
	// TODO Implement to adjust invalid solution
}

#endif /* PERMUTATIONDECODER_H_ */
