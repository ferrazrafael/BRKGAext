/*
 * main_irace.cpp
 *
 *  Created on: Dec 11, 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>


#include <iostream>
#include <iomanip>
#include "PermutationDecoder.h"
#include "MTRand.h"
#include "BRKGAext.h"

using namespace std;

int main(int argc, char* argv[]) {
	unsigned n, p;
	double pe, pm, rhoe;
	unsigned K, X_INTVL, X_NUMBER, MAX_GENS;
	const unsigned MAXT = 2;	// number of threads for parallel decoding

	vector<string> arguments(argv + 1, argv + argc);
	for(unsigned i=1; i<arguments.size(); i+=2)
	{
		string::size_type sz;
		if(arguments[i]== "--p")
			p = n*stoi(arguments[i+1], &sz);
		else if(arguments[i]== "--pe")
			pe = stod(arguments[i+1], &sz);
		else if(arguments[i]== "--pm")
			pm = stod(arguments[i+1], &sz);
		else if(arguments[i]== "--rhoe")
			rhoe = stod(arguments[i+1], &sz);
		else if(arguments[i]== "--K")
			K = stoi(arguments[i+1], &sz);
		else if(arguments[i]== "--X_INTVL")
			X_INTVL = stoi(arguments[i+1], &sz);
		else if(arguments[i]== "--MAX_GENS")
			MAX_GENS = stoi(arguments[i+1], &sz);
		else if(arguments[i]== "--X_NUMBER")
			X_NUMBER = stoi(arguments[i+1], &sz);
	}

	const long unsigned rngSeed = 0;	// seed to the random number generator
	MTRand rng(rngSeed);				// initialize the random number generator

	vector<unsigned> items = {1, 2, 3, 4, 5};     // items of the permutation
	Permutation permutation(items, items.size());
	PermutationDecoder decoder(permutation);

	BRKGAext< PermutationDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT); // initialize the decoder

	//algorithm.useDynamicOperators(true);
	//algorithm.useEliteDiversification(true, 0.1);
	//knapsack.initializeNonRandom(func);

	unsigned generation = 0;		// current generation
	do {
		algorithm.evolve();	// evolve the population for one generation

		if((++generation) % X_INTVL == 0) {
			algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
		}

		//algorithm.applyHeuristic<unsigned>( bind(&PermutationDecoder::twoSwap, decoder, placeholders::_1) );
	} while (generation < MAX_GENS);

	return unsigned( algorithm.getBestFitness() );
}*/



