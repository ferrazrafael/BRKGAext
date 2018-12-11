/*
 * BRKGAext.h
 *
 *  Created on: Jun 29, 2018
 *      Author: Rafael Louback Ferraz <ferrazrafael@yahoo.com.br>
 *
 */

#ifndef BRKGAEXT_H_
#define BRKGAEXT_H_

#include <functional>
#include "BRKGA.h"

template< class Decoder, class RNG >
class BRKGAext : public BRKGA<Decoder, RNG > {
public:
	BRKGAext(unsigned n, unsigned p, double pe, double pm, double rhoe,
			const Decoder& refDecoder, RNG& refRNG, unsigned K = 1, unsigned MAX_THREADS = 1);
	virtual ~BRKGAext();

	void evolve(unsigned generations = 1);

	// DynamicOperator
	void useDynamicOperators(bool on);

	// Elite diversification
	void useEliteDiversification(bool on, double diversityRate = 0.05);

	// Apply Heuristic implemented by user
	template <typename ItemType>
	void applyHeuristic(std::function< std::vector<ItemType> (const std::vector<ItemType>&) > heuristic);

	// Initialize population based on initial solution
	void heuristicInitialization(std::function< std::vector<unsigned> (unsigned) > heuristicInit);

protected:
	// Just to shorten base class variable use
	using BRKGA<Decoder, RNG>::n;
	using BRKGA<Decoder, RNG>::p;
	using BRKGA<Decoder, RNG>::pe;
	using BRKGA<Decoder, RNG>::pm;
	using BRKGA<Decoder, RNG>::rhoe;
	using BRKGA<Decoder, RNG>::refRNG;
	using BRKGA<Decoder, RNG>::refDecoder;
	using BRKGA<Decoder, RNG>::K;
	using BRKGA<Decoder, RNG>::MAX_THREADS;
	using BRKGA<Decoder, RNG>::previous;
	using BRKGA<Decoder, RNG>::current;

	// Elite diversification
	bool eliteDiversification = false;
	double diversityRate = 0.0;

	// Adaptive parameters
	bool dynamicOperators = false;
	std::vector<unsigned> dynPe;
	std::vector<unsigned> dynPm;
	std::vector<double> dynRhoe;

	std::vector<double> generateChromosome(const std::vector<unsigned>& permutation);
	void diversifyElite();
	double eliteAverage(const Population& population, unsigned auxPe);
	double eliteStandardDeviation(const Population& population, double avg, unsigned auxPe);
	void updateDynamicOperators();
	void dynamicEvolution(Population& curr, Population& next, const unsigned k);
};

template< class Decoder, class RNG >
BRKGAext<Decoder, RNG>::BRKGAext(unsigned _n, unsigned _p, double _pe, double _pm, double _rhoe,
		const Decoder& decoder, RNG& rng, unsigned _K, unsigned MAX) :
		BRKGA<Decoder, RNG>(_n, _p, _pe, _pm, _rhoe, decoder, rng, _K, MAX) {
	dynPe.assign(K, pe); // assigns 'pe' value K times
	dynPm.assign(K, pm); // assigns 'pm' value K times
	dynRhoe.assign(K, rhoe); // assigns 'rhoe' value K times
}

template< class Decoder, class RNG >
BRKGAext<Decoder, RNG>::~BRKGAext() {
	BRKGA<Decoder, RNG>::~BRKGA();
}

template< class Decoder, class RNG >
void BRKGAext< Decoder, RNG >::evolve(unsigned generations) {
	if(dynamicOperators){
		if(generations == 0) { throw std::range_error("Cannot evolve for 0 generations."); }

		for(unsigned i = 0; i < generations; ++i) {
			for(unsigned j = 0; j < K; ++j) {
				dynamicEvolution(*(current[j]), *(previous[j]), j);	// Evolve population
				std::swap(current[j], previous[j]);		// Update generation
			}
		}
		updateDynamicOperators();
	}
	else {
		BRKGA<Decoder, RNG>::evolve(generations);
	}

	if(eliteDiversification) {
		diversifyElite();
	}
}

template<class Decoder, class RNG>
void BRKGAext<Decoder, RNG>::useDynamicOperators(bool on) {
	dynamicOperators = on;
}

template<class Decoder, class RNG>
void BRKGAext<Decoder, RNG>::useEliteDiversification(bool on, double _diversityRate) {
	eliteDiversification = on;
	diversityRate = _diversityRate;
}

template<class Decoder, class RNG>
template<typename ItemType>
void BRKGAext<Decoder, RNG>::applyHeuristic(std::function< std::vector<ItemType> (const std::vector<ItemType>&) > heuristic) {
	for(int i = 0; i < int(K); ++i) {
		for(int j = 0; j < int(p); ++j) {
			std::vector<ItemType> solution = refDecoder.decodeSolution((*current[i])(j)) ;
			solution = heuristic(solution);
			refDecoder.correctChromosome((*current[i])(j), solution);
			current[i]->setFitness(j, refDecoder.decode((*current[i])(j)));
		}
	}
}

template< class Decoder, class RNG >
void BRKGAext< Decoder, RNG >::heuristicInitialization(std::function<std::vector<unsigned> (unsigned)> heuristicInit) {
	for(unsigned i = 0; i < K; ++i) {
		for(unsigned j = 0; j < p; ++j) {
			std::vector<unsigned> permutation = heuristicInit(n);
			if(permutation.size() != n) {
				throw std::range_error("Permutation size differs from chromosome size.");
			}

			(*current[i])(j) = generateChromosome(permutation);
		}

		// Decode:
#ifdef _OPENMP
#pragma omp parallel for num_threads(MAX_THREADS)
#endif
		for(int j = 0; j < int(p); ++j) {
			current[i]->setFitness(j, refDecoder.decode((*current[i])(j)) );
		}

		// Sort:
		current[i]->sortFitness();
	}
}

template< class Decoder, class RNG >
std::vector<double> BRKGAext< Decoder, RNG >::generateChromosome(const std::vector<unsigned>& permutation) {
	std::vector<double> chromosome(n);
	std::vector<double> keys(n);

	for(unsigned i = 0; i < n; ++i) {
		keys[i] = refRNG.rand();
	}
	std::sort(keys.begin(), keys.end());

	for(unsigned i = 0; i < n; ++i) {
		chromosome[i] = keys[permutation[i]];
	}

	return chromosome;
}

template< class Decoder, class RNG >
void BRKGAext<Decoder, RNG>::diversifyElite() {
	// for each population check if elite is to homogeneous
	for(unsigned i = 0; i < K; ++i) {
		// auxiliary 'pe' that uses dynamic pe if dynamicOperators are on or default pe otherwise
		unsigned _pe = dynamicOperators ? dynPe[i] : pe;

		double avg = eliteAverage(*current[i], _pe);
		double percent = (eliteStandardDeviation(*current[i], avg, _pe) * 100.0) / current[i]->getBestFitness() / 100.0;

		if(percent < diversityRate) {
			// replace elite chromosomes that are below diversity rate with mutants
			for(unsigned j = 1; j < _pe; ++j) {
				double percentualDeviation = fabs( (current[i]->getFitness(j) - avg) / avg );

				if(percentualDeviation < diversityRate) {
					// generating mutant chromosome
					for(unsigned k = 0; k < n; ++k) { (*current[i])(j, k) = refRNG.rand(); }
				}
			}
			// Time to compute fitness, in parallel:
#ifdef _OPENMP
#pragma omp parallel for num_threads(MAX_THREADS)
#endif
			for(int j = 1; j < int(_pe); ++j) {
				current[i]->setFitness( j, refDecoder.decode(current[i]->population[i]) );
			}

			// Now we must sort 'current' by fitness, since things might have changed:
			current[i]->sortFitness();
		}
	}
}

template< class Decoder, class RNG >
double BRKGAext<Decoder, RNG>::eliteAverage(const Population& population, unsigned _pe) {
	double avg = 0.0;
	for(unsigned i = 0; i < _pe; ++i) {
		avg += population.getFitness(i);
	}
	avg /= double(_pe); // divide sum of elite fitness by auxiliary 'pe'
	return avg;
}

template< class Decoder, class RNG >
double BRKGAext<Decoder, RNG>::eliteStandardDeviation(const Population& population, double avg, unsigned _pe) {
	double sum = 0.0;
	for(unsigned i = 0; i < _pe; ++i) {
		sum += population.getFitness(i) - avg;
	}
	sum /= double(_pe); // divide sum of elite fitness by auxiliary 'pe'

	return sqrt(sum);
}

template< class Decoder, class RNG >
void BRKGAext< Decoder, RNG >::dynamicEvolution(Population& curr, Population& next, const unsigned k) {
	// We now will set every chromosome of 'current', iterating with 'i':
	unsigned i = 0;	// Iterate chromosome by chromosome
	unsigned j = 0;	// Iterate allele by allele

	// 2. The 'pe' best chromosomes are maintained, so we just copy these into 'current':
	while(i < dynPe[k]) {
		for(j = 0 ; j < n; ++j) { next(i,j) = curr(curr.fitness[i].second, j); }

		next.fitness[i].first = curr.fitness[i].first;
		next.fitness[i].second = i;
		++i;
	}

	// 3. We'll mate 'p - pe - pm' pairs; initially, i = pe, so we need to iterate until i < p - pm:
	while(i < p - dynPm[k]) {
		// Select an elite parent:
		const unsigned eliteParent = (refRNG.randInt(dynPe[k] - 1));

		// Select a non-elite parent:
		const unsigned noneliteParent = dynPe[k] + (refRNG.randInt(p - dynPe[k] - 1));

		// Mate:
		for(j = 0; j < n; ++j) {
			const unsigned sourceParent = ((refRNG.rand() < dynRhoe[k]) ? eliteParent : noneliteParent);

			next(i, j) = curr(curr.fitness[sourceParent].second, j);
		}

		++i;
	}

	// We'll introduce 'pm' mutants:
	while(i < p) {
		for(j = 0; j < n; ++j) { next(i, j) = refRNG.rand(); }
		++i;
	}

	// Time to compute fitness, in parallel:
#ifdef _OPENMP
#pragma omp parallel for num_threads(MAX_THREADS)
#endif
	for(int i = int(dynPe[k]); i < int(p); ++i) {
		next.setFitness( i, refDecoder.decode(next.population[i]) );
	}

	// Now we must sort 'current' by fitness, since things might have changed:
	next.sortFitness();
}

template< class Decoder, class RNG >
void BRKGAext<Decoder, RNG>::updateDynamicOperators() {
	// update operators
	for(unsigned i = 0; i < K; ++i) {
		double avgFitness = 0.0;
		for(unsigned j = 0; j < p; ++j) {
			avgFitness += current[i]->getFitness(j);
		}
		avgFitness /= double(p); // divide sum of all fitness by p

		double bestFitness = current[i]->getBestFitness();

		// generating adaptive pe's
		double k = double(pe) / 2;
		double avg = 0.0;
		for(unsigned j = 0; j < p; ++j) {
			avg += fabs( k * (current[i]->getFitness(j) - bestFitness) / (avgFitness - bestFitness) );
		}
		avg /= double(p); // divide sum of all Pm's by p
		dynPe[i] = unsigned(avg);

		// generating adaptive pm's
		k = double(pm);
		avg = 0.0;
		for(unsigned j = 0; j < p; ++j) {
			avg += fabs( k / (current[i]->getFitness(j) - bestFitness) ) ;
		}
		avg /= double(p); // divide sum of all Pm's by p
		dynPm[i] = unsigned(avg);

		// generating adaptive rhoe's
		avg = 0.0;
		for(unsigned j = 0; j < p; ++j) {
			avg += 0.5 + fabs( 0.25 * (current[i]->getFitness(j) - bestFitness) / (avgFitness - bestFitness) );
		}
		dynRhoe[i] = avg / double(p); // divide sum of all Rhoe's by p
	}
}

#endif /* BRKGAEXT_H_ */
