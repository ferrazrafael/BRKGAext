#include <iostream>
#include <iomanip>
#include "KnapsackDecoder.h"
#include "PermutationDecoder.h"
#include "MTRand.h"
#include "BRKGAext.h"

using namespace std;

int main(int argc, char* argv[]) {
	const unsigned n = 4;		// size of chromosomes
	const unsigned p = 1000;	// size of population
	const double pe = 0.20;		// fraction of population to be the elite-set
	const double pm = 0.10;		// fraction of population to be replaced by mutants
	const double rhoe = 0.70;	// probability that offspring inherit an allele from elite parent
	const unsigned K = 3;		// number of independent populations
	const unsigned MAXT = 2;	// number of threads for parallel decoding

	const long unsigned rngSeed = 0;	// seed to the random number generator
	MTRand rng(rngSeed);				// initialize the random number generator

	vector<double> values = {30, 14, 16, 9};      // value of each item on the knapsack
	vector<double> weights = {6, 3, 4, 2};     // weight of each item on the knapsack
	double W = 100;
	Knapsack knapsack(n, values, weights, W);
	KnapsackDecoder decoder(knapsack);

	//vector<unsigned> items = {1, 2, 3, 4, 5};     // items of the permutation
	//Permutation permutation(items, items.size());
	//PermutationDecoder decoder(permutation);


	BRKGAext< KnapsackDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT); // initialize the decoder

	algorithm.useDynamicOperators(true);
	//algorithm.useEliteDiversification(0.1);
	//knapsack.initializeNonRandom(func);

	unsigned generation = 0;		// current generation
	const unsigned X_INTVL = 100;	// exchange best individuals at every 100 generations
	const unsigned X_NUMBER = 2;	// exchange top 2 best
	const unsigned MAX_GENS = 1000;	// run for 1000 gens
	do {
		algorithm.evolve();	// evolve the population for one generation

		if((++generation) % X_INTVL == 0) {
			algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
		}

		//algorithm.applyHeuristic<unsigned>( bind(&PermutationDecoder::twoSwap, decoder, placeholders::_1) );
	} while (generation < MAX_GENS);

	cout << "Best solution found has objective value = "
			<< fixed << std::setw(11)
	<< setprecision(6)
	<< algorithm.getBestFitness() << endl;


	return 0;
}

std::vector<unsigned> func(unsigned n){
	std::vector<unsigned> permutation(n);
	for(unsigned i = n; i > 0; --i)
		permutation[i-1] = n-i;

	return permutation;
};

