/*
 * Knapsack.h
 *
 *  Created on: Oct 4, 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>
 */

#ifndef KNAPSACK_H_
#define KNAPSACK_H_

#include <vector>

class Knapsack {
public:
	Knapsack(unsigned n, const std::vector<double>& values, const std::vector<double>& weights, double W);
	virtual ~Knapsack();

	unsigned getN() const;
	const std::vector<double>& getValues() const;
	double getValue(unsigned i) const;
	double getW() const;
	const std::vector<double>& getWeights() const;
	double getWeight(unsigned i) const;



protected:
	unsigned n; // number of items
	std::vector<double> values; // value of each item 'i'
	std::vector<double> weights; // weight of each item 'i'
	double W; // max weight supported by the Knapsack
};

// Local Search
std::vector<bool> twoSwap(const std::vector<bool>& selection, const Knapsack& knapsack);

#endif /* KNAPSACK_H_ */
