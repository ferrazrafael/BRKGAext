/*
 * PermutationDecoder.h
 *
 *  Created on: Nov 20, 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>
 */

#ifndef PERMUTATIONDECODER_H_
#define PERMUTATIONDECODER_H_

#include <vector>
#include "Permutation.h"

class PermutationDecoder {
public:
	PermutationDecoder(const Permutation permutation);
	virtual ~PermutationDecoder();

	double decode(std::vector<double>& chromosome) const;
	double computeFitness(const std::vector<ItemType>& solution) const;
	std::vector<ItemType> decodeSolution(std::vector<double>& chromosome) const;
	void correctChromosome(std::vector<double>& chromosome, const std::vector<ItemType>& solution) const;
	bool isValid(const std::vector<ItemType>& solution);

	// Local Search
	std::vector<ItemType> twoSwap(const std::vector<ItemType>& solution);
	std::vector<ItemType> kSwap(const std::vector<ItemType>& solution, unsigned k);
	std::vector<ItemType> bestInsertion(const std::vector<ItemType>& solution);
	std::vector<ItemType> twoOpt(const std::vector<ItemType>& solution);

protected:
	const Permutation permutation;
	void adjustSolution(std::vector<ItemType>& solution) const;
};

#endif /* PERMUTATIONDECODER_H_ */
