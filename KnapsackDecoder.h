/*
 * Knapsack.h
 *
 *  Created on: 26 de set de 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>
 */

#ifndef KNAPSACKDECODER_H_
#define KNAPSACKDECODER_H_

#include <vector>

class KnapsackDecoder {
public:
	KnapsackDecoder(unsigned n, std::vector<double> values, std::vector<double> weights, double W);
	virtual ~KnapsackDecoder();

	double decode(std::vector<double>& chromosome) const;
	unsigned getN() const;
	const std::vector<double>& getValues() const;
	double getW() const;
	const std::vector<double>& getWeights() const;

protected:
	unsigned n; // number of items
	std::vector<double> values; // value of each item 'i'
	std::vector<double> weights; // weight of each item 'i'
	double W; // max weight supported by the Knapsack

	void adjustSelection(std::vector<bool>& selection) const;
	void correctChromosome(std::vector<double>& chromosome, std::vector<bool>& selection) const;
};

#endif /* KNAPSACKDECODER_H_ */
