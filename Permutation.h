/*
 * Permutation.h
 *
 *  Created on: Nov 20, 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>
 */

#ifndef PERMUTATION_H_
#define PERMUTATION_H_

#include <vector>

typedef unsigned ItemType;

class Permutation {
public:
	Permutation(std::vector<ItemType> items, unsigned n) : items(items), n(n) {};
	virtual ~Permutation() {};

	const std::vector<ItemType>& getItems() const { return items; };
	ItemType getItem(unsigned i) const { return items[i]; };
	unsigned getN() const { return n; };

protected:
	std::vector<ItemType> items;
	unsigned n;
};

#endif /* PERMUTATION_H_ */
