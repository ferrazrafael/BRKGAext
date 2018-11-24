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
	double computeFitness(const std::vector<unsigned>& solution) const;
	std::vector<unsigned> decodeSolution(std::vector<double>& chromosome) const;
	void correctChromosome(std::vector<double>& chromosome, std::vector<unsigned>& solution) const;

	// Local Search
	std::vector<unsigned> twoSwap(const std::vector<unsigned>& solution);

protected:
	const Permutation permutation;
	void adjustSolution(std::vector<unsigned>& solution) const;
};

#endif /* PERMUTATIONDECODER_H_ */
