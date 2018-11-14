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

	// Adaptive parameters
	void useAdaptiveParameters();

	// Elite diversification
	void useEliteDiversification(double diversityRate = 0.05);

	// Initialize population based on initial solution
	void initializeNonRandom(std::function<std::vector<unsigned> (unsigned)> initFunc);

	// Apply Heuristic implemented by user
	template <typename SolutionType>
	void applyHeuristic(std::function< SolutionType (const SolutionType&) > heuristic);

protected:
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
	bool adaptiveParameters = false;
	std::vector<unsigned> adpPe;
	std::vector<unsigned> adpPm;
	std::vector<double> adpRhoe;

	void diversifyElite();
	void updateAdaptiveParameters();
	void adaptiveEvolution(Population& curr, Population& next, const unsigned k);
	double average(const Population& population);
	double standardDeviation(const Population& population, double avg);
	std::vector<double> generateChromosome(const std::vector<unsigned>& permutation);
};

template< class Decoder, class RNG >
BRKGAext<Decoder, RNG>::BRKGAext(unsigned _n, unsigned _p, double _pe, double _pm, double _rhoe,
		const Decoder& decoder, RNG& rng, unsigned _K, unsigned MAX) :
		BRKGA<Decoder, RNG>(_n, _p, _pe, _pm, _rhoe, decoder, rng, _K, MAX) {
	adpPe.assign(K, pe);
	adpPm.assign(K, pm);
	adpRhoe.assign(K, rhoe);
}

template< class Decoder, class RNG >
BRKGAext<Decoder, RNG>::~BRKGAext() {
	BRKGA<Decoder, RNG>::~BRKGA();
}

template< class Decoder, class RNG >
void BRKGAext< Decoder, RNG >::evolve(unsigned generations) {
	if(adaptiveParameters){
		if(generations == 0) { throw std::range_error("Cannot evolve for 0 generations."); }

		for(unsigned i = 0; i < generations; ++i) {
			for(unsigned j = 0; j < K; ++j) {
				adaptiveEvolution(*(current[j]), *(previous[j]), j);	// First evolve the population (curr, next)
				std::swap(current[j], previous[j]);		// Update (prev = curr; curr = prev == next)
			}
		}
		updateAdaptiveParameters();
	}
	else {
		BRKGA<Decoder, RNG>::evolve(generations);
	}

	if(eliteDiversification) {
		diversifyElite();
	}
}

template< class Decoder, class RNG >
inline void BRKGAext< Decoder, RNG >::adaptiveEvolution(Population& curr, Population& next, const unsigned k) {
	// We now will set every chromosome of 'current', iterating with 'i':
	unsigned i = 0;	// Iterate chromosome by chromosome
	unsigned j = 0;	// Iterate allele by allele

	// 2. The 'pe' best chromosomes are maintained, so we just copy these into 'current':
	while(i < adpPe[k]) {
		for(j = 0 ; j < n; ++j) { next(i,j) = curr(curr.fitness[i].second, j); }

		next.fitness[i].first = curr.fitness[i].first;
		next.fitness[i].second = i;
		++i;
	}

	// 3. We'll mate 'p - pe - pm' pairs; initially, i = pe, so we need to iterate until i < p - pm:
	while(i < p - adpPm[k]) {
		// Select an elite parent:
		const unsigned eliteParent = (refRNG.randInt(adpPe[k] - 1));

		// Select a non-elite parent:
		const unsigned noneliteParent = adpPe[k] + (refRNG.randInt(p - adpPe[k] - 1));

		// Mate:
		for(j = 0; j < n; ++j) {
			const unsigned sourceParent = ((refRNG.rand() < adpRhoe[k]) ? eliteParent : noneliteParent);

			next(i, j) = curr(curr.fitness[sourceParent].second, j);

			//next(i, j) = (refRNG.rand() < rhoe) ? curr(curr.fitness[eliteParent].second, j) :
			//		                              curr(curr.fitness[noneliteParent].second, j);
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
	for(int i = int(adpPe[k]); i < int(p); ++i) {
		next.setFitness( i, refDecoder.decode(next.population[i]) );
	}

	// Now we must sort 'current' by fitness, since things might have changed:
	next.sortFitness();
}

template<class Decoder, class RNG>
void BRKGAext<Decoder, RNG>::useAdaptiveParameters() {
	adaptiveParameters = true;
}

template<class Decoder, class RNG>
void BRKGAext<Decoder, RNG>::useEliteDiversification(double _diversityRate) {
	eliteDiversification = true;
	diversityRate = _diversityRate;
}

template< class Decoder, class RNG >
double BRKGAext<Decoder, RNG>::average(const Population& population) {
	double avg = 0.0;
	for(unsigned i = 0; i < pe; ++i) {
		avg += population.getFitness(i);
	}
	avg /= double(pe); // divide sum of all fitness by p
	return avg;
}

template< class Decoder, class RNG >
double BRKGAext<Decoder, RNG>::standardDeviation(const Population& population, double avg) {
	double sum = 0.0;
	for(unsigned i = 0; i < pe; ++i) {
		sum += population.getFitness(i) - avg;
	}
	sum /= double(pe); // divide sum of all fitness by p

	return sqrt(sum);
}

template< class Decoder, class RNG >
void BRKGAext<Decoder, RNG>::diversifyElite() {
	// for each population check if elite is to homogeneous
	for(unsigned i = 0; i < K; ++i) {
		double avg = average(*current[i]);
		double percent = (standardDeviation(*current[i], avg) * 100.0) / current[i]->getBestFitness() / 100.0;

		if(percent < diversityRate) {
			// replace elite chromosomes that are below diversity rate with mutants
			for(unsigned j = 1; j < pe; ++j) {
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
			for(int j = 1; j < int(pe); ++j) {
				current[i]->setFitness( j, refDecoder.decode(current[i]->population[i]) );
			}

			// Now we must sort 'current' by fitness, since things might have changed:
			current[i]->sortFitness();
		}
	}
}

template< class Decoder, class RNG >
void BRKGAext<Decoder, RNG>::updateAdaptiveParameters() {
	// update parameters
	for(unsigned i = 0; i < K; ++i) {
		double avgFitness = 0.0;
		for(unsigned j = 0; j < p; ++j) {
			avgFitness += current[i]->getFitness(j);
		}
		avgFitness /= double(p); // divide sum of all fitness by p

		double bestFitness = current[i]->getBestFitness();

		// generating adaptive pm's
		double k = (pm / double(p)) / 2;
		double avg = 0.0;
		for(unsigned j = 0; j < p; ++j) {
			avg += fabs(k * (current[i]->getFitness(j) - bestFitness) / (avgFitness - bestFitness)) ;
		}
		avg /= double(p); // divide sum of all Pm's by p
		adpPm[i] = unsigned(p * avg);

		// generating adaptive pe's
		k = (pe / double(p)) / 2;
		avg = 0.0;
		for(unsigned j = 0; j < p; ++j) {
			avg += fabs( k * (current[i]->getFitness(j) - bestFitness) / (avgFitness - bestFitness) );
		}
		avg /= double(p); // divide sum of all Pm's by p
		adpPe[i] = unsigned(p * avg);

		// generating adaptive rhoe's
		avg = 0.0;
		for(unsigned j = 0; j < p; ++j) {
			avg += 0.5 + fabs( 0.25 * (current[i]->getFitness(j) - bestFitness) / (avgFitness - bestFitness) );
		}
		adpRhoe[i] = avg / double(p); // divide sum of all Rhoe's by p

		/*std::cout << "Pm: " << adpPm[i] << "  ";
		std::cout << "Pe: " << adpPe[i] << "  ";
		std::cout << "Rhoe: " << adpRhoe[i] << std::endl;*/
	}
}

template<class Decoder, class RNG>
template<typename SolutionType>
void BRKGAext<Decoder, RNG>::applyHeuristic(std::function< SolutionType (const SolutionType&) > heuristic) {
	for(int i = 0; i < int(K); ++i) {
		for(int j = 0; j < int(p); ++j) {
			SolutionType solution = refDecoder.computeSolution((*current[i])(j)) ;
			solution = heuristic(solution);
			refDecoder.correctChromosome((*current[i])(j), solution);
			current[i]->setFitness(j, refDecoder.decode((*current[i])(j)));
		}
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
void BRKGAext< Decoder, RNG >::initializeNonRandom(std::function<std::vector<unsigned> (unsigned)> initFunc) {
	using std::range_error;

	for(unsigned i = 0; i < K; ++i) {
		for(unsigned j = 0; j < p; ++j) {
			std::vector<unsigned> permutation = initFunc(n);
			if(permutation.size() != n) {
				throw range_error("Permutation size differs from chromosome size.");
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

#endif /* BRKGAEXT_H_ */
