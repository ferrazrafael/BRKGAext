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
	Permutation();
	virtual ~Permutation();
	const std::vector<ItemType>& getItems() const;
	ItemType getItem(unsigned i) const;
	unsigned getN() const;

protected:
	std::vector<ItemType> items;
	unsigned n;
};

#endif /* PERMUTATION_H_ */
