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
	vector<unsigned> solution  = decodeSolution(chromosome);

	return computeFitness(solution);
}

double PermutationDecoder::computeFitness(const std::vector<unsigned>& solution) const {
	// TODO API user implements fitness computation based on problem data
}

vector<unsigned> PermutationDecoder::decodeSolution(vector<double>& chromosome) const {
}

void PermutationDecoder::correctChromosome(vector<double>& chromosome, vector<unsigned>& solution) const {
}

vector<unsigned> PermutationDecoder::twoSwap(const vector<unsigned>& solution) {
}

void PermutationDecoder::adjustSolution(vector<unsigned>& solution) const {
}
