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
	Knapsack(unsigned n, const std::vector<double>& values, const std::vector<double>& weights, double W)
			: n(n), values(values), weights(weights), W(W) {};
	virtual ~Knapsack() {};

	unsigned getN() const { return n; };
	const std::vector<double>& getValues() const { return values; };
	double getValue(unsigned i) const { return values[i]; };
	double getW() const { return W; };
	const std::vector<double>& getWeights() const { return weights; };
	double getWeight(unsigned i) const { return weights[i]; };

protected:
	unsigned n; // number of items
	std::vector<double> values; // value of each item 'i'
	std::vector<double> weights; // weight of each item 'i'
	double W; // max weight supported by the Knapsack
};



#endif /* KNAPSACK_H_ */
