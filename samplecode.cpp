#include <iostream>
#include <iomanip>
#include "KnapsackDecoder.h"
#include "MTRand.h"
#include "BRKGAext.h"

using namespace std;

int main(int argc, char* argv[]) {
	const unsigned n = 100;		// size of chromosomes
	const unsigned p = 1000;	// size of population
	const double pe = 0.20;		// fraction of population to be the elite-set
	const double pm = 0.10;		// fraction of population to be replaced by mutants
	const double rhoe = 0.70;	// probability that offspring inherit an allele from elite parent
	const unsigned K = 3;		// number of independent populations
	const unsigned MAXT = 2;	// number of threads for parallel decoding

	vector<double> values = {30, 14, 16, 9};      // value of each item on the knapsack
	vector<double> weights = {6, 3, 4, 2};     // weight of each item on the knapsack
	double W = 1;

	const long unsigned rngSeed = 0;	// seed to the random number generator
	MTRand rng(rngSeed);				// initialize the random number generator

	Knapsack knapsack(n, values, weights, W);
	KnapsackDecoder decoder(knapsack);

	BRKGAext< KnapsackDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT); // initialize the decoder

	algorithm.useAdaptiveParameters();
	algorithm.useEliteDiversification(0.1);
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

		algorithm.applyHeuristic<vector<bool>, Knapsack>(twoSwap, knapsack);
	} while (generation < MAX_GENS);

	cout << "Best solution found has objective value = "
	     << fixed << std::setw(11)
	     << setprecision(6)
	     << algorithm.getBestFitness() << endl;


	return 0;
}

/*std::vector<unsigned> func(unsigned n){
	std::vector<unsigned> permutation(n);
	for(unsigned i = n; i > 0; --i)
		permutation[i-1] = n-i;

	return permutation;
};

int main(int argc, char* argv[]) {
	const unsigned n = 100;		// size of chromosomes
	const unsigned p = 1000;	// size of population
	const double pe = 0.20;		// fraction of population to be the elite-set
	const double pm = 0.10;		// fraction of population to be replaced by mutants
	const double rhoe = 0.70;	// probability that offspring inherit an allele from elite parent
	const unsigned K = 3;		// number of independent populations
	const unsigned MAXT = 2;	// number of threads for parallel decoding

	SampleDecoder decoder;			// initialize the decoder

	const long unsigned rngSeed = 0;	// seed to the random number generator
	MTRand rng(rngSeed);				// initialize the random number generator

	// initialize the BRKGA-based heuristic
	BRKGAext< SampleDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);
	algorithm.useAdaptiveParameters();
	algorithm.useEliteDiversification(0.1);
	//algorithm.initializeNonRandom(func);

	unsigned generation = 0;		// current generation
	const unsigned X_INTVL = 100;	// exchange best individuals at every 100 generations
	const unsigned X_NUMBER = 2;	// exchange top 2 best
	const unsigned MAX_GENS = 1000;	// run for 1000 gens
	do {
		algorithm.evolve();	// evolve the population for one generation

		if((++generation) % X_INTVL == 0) {
			algorithm.exchangeElite(X_NUMBER);	// exchange top individuals
		}
	} while (generation < MAX_GENS);

	std::cout << "Best solution found has objective value = "
	 		  << std::fixed << std::setw(11)
    		  << std::setprecision(6)
	          << algorithm.getBestFitness() << std::endl;

	return 0;
}*/
