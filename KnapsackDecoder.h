/*
 * Knapsack.h
 *
 *  Created on: Sep 26, 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>
 */

#ifndef KNAPSACKDECODER_H_
#define KNAPSACKDECODER_H_

#include <vector>
#include "Knapsack.h"

class KnapsackDecoder {
public:
	KnapsackDecoder(const Knapsack& knapsack);
	virtual ~KnapsackDecoder();

	double decode(std::vector<double>& chromosome) const;
	double computeFitness(const std::vector<bool>& solution) const;
	std::vector<bool> decodeSolution(std::vector<double>& chromosome) const;
	void correctChromosome(std::vector<double>& chromosome, std::vector<bool>& solution) const;

	// Local Search
	std::vector<bool> twoSwap(const std::vector<bool>& selection);

protected:
	const Knapsack knapsack;
	void adjustSolution(std::vector<bool>& solution) const;

};


#endif /* KNAPSACKDECODER_H_ */
