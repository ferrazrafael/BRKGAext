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
	std::vector<bool> getSolution(std::vector<double>& chromosome) const;
	void correctChromosome(std::vector<double>& chromosome, std::vector<bool>& selection) const;

protected:
	const Knapsack knapsack;
	void adjustSelection(std::vector<bool>& selection) const;

};


#endif /* KNAPSACKDECODER_H_ */
