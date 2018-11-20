/*
 * PermutationDecoder.h
 *
 *  Created on: 20 de nov de 2018
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
	double computeFitness(const std::vector<bool>& solution) const;
	std::vector<bool> decodeSolution(std::vector<double>& chromosome) const;
	void correctChromosome(std::vector<double>& chromosome, std::vector<bool>& solution) const;

	// Local Search
	std::vector<bool> twoSwap(const std::vector<bool>& selection);

protected:
	const Permutation permutation;
	void adjustSolution(std::vector<bool>& solution) const;
};

#endif /* PERMUTATIONDECODER_H_ */
