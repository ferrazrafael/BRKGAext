/*
 * Permutation.cpp
 *
 *  Created on: Nov 20, 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>
 */

#include "Permutation.h"

Permutation::Permutation() {
	// TODO Auto-generated constructor stub

}

Permutation::~Permutation() {
	// TODO Auto-generated destructor stub
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
