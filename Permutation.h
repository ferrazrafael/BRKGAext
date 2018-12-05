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
	Permutation(std::vector<ItemType> items, unsigned n);
	virtual ~Permutation();

	const std::vector<ItemType>& getItems() const;
	ItemType getItem(unsigned i) const;
	unsigned getN() const;

protected:
	std::vector<ItemType> items;
	unsigned n;
};

Permutation::Permutation(std::vector<ItemType> _items, unsigned _n) : items(_items), n(_n) {
}

Permutation::~Permutation() {
}

const std::vector<ItemType>& Permutation::getItems() const {
	return items;
}

ItemType Permutation::getItem(unsigned i) const {
	return items[i];
}

unsigned Permutation::getN() const {
	return n;
}

#endif /* PERMUTATION_H_ */
