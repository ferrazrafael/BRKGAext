/*
 * PermutationDecoder.cpp
 *
 *  Created on: Nov 20, 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>
 */

#include "PermutationDecoder.h"

using namespace std;

PermutationDecoder::PermutationDecoder(const Permutation _permutation) : permutation(_permutation) {
}

PermutationDecoder::~PermutationDecoder() {
}

double PermutationDecoder::decode(std::vector<double>& chromosome) const {
}

double PermutationDecoder::computeFitness(const std::vector<bool>& solution) const {
	// TODO API user implements fitness computation based on problem data
}

vector<bool> PermutationDecoder::decodeSolution(vector<double>& chromosome) const {
}

void PermutationDecoder::correctChromosome(vector<double>& chromosome, vector<bool>& solution) const {
}

vector<bool> PermutationDecoder::twoSwap(const vector<bool>& selection) {
}

void PermutationDecoder::adjustSolution(vector<bool>& solution) const {
}
