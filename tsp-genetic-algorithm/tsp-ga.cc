#include "tsp-ga.hh"
#include <iostream>
#include <algorithm> 	// std::random_shuffle(), std::iter_swap(), std::sort()
#include <set>			// std::set
#include <cstdlib>  	// rand()

using std::vector;
using std::set;
using std::cout;
using std::endl;

/**
 * @description: generate an random numPoints length's Genome by random shuffle.
 * @param {int} numPoints: Genome's gene number.
 * @return {*}
 */
TSPGenome::TSPGenome(const int numPoints) {
	for(int i = 0; i < numPoints; ++i)
		tspgGenome.push_back(i);

	std::random_shuffle (tspgGenome.begin(), tspgGenome.end());
	tspgCircuitLength = -1;
}

/**
 * @description: Initializes the genome from the specified visit order
 * @param {vector<int>} &order 
 * @return {*}
 */
TSPGenome::TSPGenome(const vector<int> &order) {
	tspgGenome = order;
	tspgCircuitLength = -1;
}

TSPGenome::~TSPGenome(){
	// no-op
}

/**
 * @description: compute the circuit length by specified order.
 * @param {vector<Point>} &points: the vector of 3-D coordinate.
 * @return {*}
 */
void TSPGenome::computeCircuitLength(const vector<Point> &points) {
	size_t sizeGenome = tspgGenome.size();
	tspgCircuitLength = points[tspgGenome.front()].distanceTo(points[tspgGenome.back()]);

	for(unsigned i = 0, j = 1; i < sizeGenome - 1; ++i, ++j) {
		tspgCircuitLength += points[tspgGenome[i]].distanceTo(points[tspgGenome[j]]);
	}
}

/**
 * @description: acquire the genome orders.
 * @return {*}
 */
vector<int> TSPGenome::getOrder() const {
	return tspgGenome;
}

double TSPGenome::getCircuitLength() const {
	return tspgCircuitLength;
}

/**
 * @description: mutates the genome by swapping two randomly selected values in the order vector.
 * @return {*}
 */
void TSPGenome::mutate() {
	unsigned swapIndexA;
	unsigned swapIndexB;

	swapIndexA = rand() % tspgGenome.size();
	do {
		swapIndexB = rand() % tspgGenome.size();
	} while (swapIndexA == swapIndexB);

	std::iter_swap(tspgGenome.begin() + swapIndexA, tspgGenome.begin() + swapIndexB);
}

/**
 * @description: used to vary the programming of a chromosome from one generation to the next.
 * @param {TSPGenome} &GA
 * @param {TSPGenome} &GB
 * @return {*}
 */
TSPGenome crossLink(const TSPGenome &GA, const TSPGenome &GB) {
	vector<int> _GA = GA.getOrder();
	vector<int> _GB = GB.getOrder();
	vector<int> offspringGenome;
	set<int> keepTrack;

	size_t sizeGenome = _GA.size();
	// generate a random index in range [2, .. , sizeGenome - 2]
	unsigned randomIndex = (rand() % (sizeGenome - 3)) + 2;

	// init offspringGenome and keepTrack
	for(unsigned i = 0; i < randomIndex; ++i) {
		offspringGenome.push_back(_GA.at(i));
		keepTrack.insert(_GA.at(i));
	}

	// crosslinking
	/**
	 * @description: 采用的交叉方式为单点交叉的方式,由于一个地方只能访问一次所以需要进行去重
	 * @return {*}
	 */
	for(const auto& gene : _GB) {
		bool found = keepTrack.find(gene) != keepTrack.end();
		if (!found) {
			offspringGenome.push_back(gene);
		}
	}

	// generate TSPGenome object from offspringGenome the return object.
	TSPGenome Offspring(offspringGenome);
	return Offspring;
}

/**
 * @description: 寻找最优路线,迭代次数为numGenerations
 * @param: numGenerations 种群迭代数
 * @param: populationSize 种群总数
 * @param: keepPopulation 交叉的个体数
 * @param: numMutations   变异的个体数
 * @return {*}
 */
TSPGenome findAShortPath(const vector<Point> &points,
                           int populationSize, int numGenerations,
                           int keepPopulation, int numMutations) {
	// Generate an initial population of random genomes
	vector<TSPGenome> population;
	size_t sizeGenome = points.size();

	for(int i = 0; i < populationSize; ++i) {
		TSPGenome genome(sizeGenome);
		population.push_back(genome);
	}

	for(int gen = 0; gen < numGenerations * 10; ++gen){

		// Compute the circuit length for each member.
		for(TSPGenome& genome : population){
			genome.computeCircuitLength(points);
		}

		/* sort the TSP genomes by fitness. 
			(The "fittest" ones should be at the lowest indexes.) */
		std::sort(population.begin(), population.end(), isShorterPath);

		/*keep the top N fittest members of the population, and replace
			the remaining members with new genomes produced from the fittest ones (crosslinking) */

		for(int i = keepPopulation; i < populationSize; ++i){

			unsigned randomIndexA;
			unsigned randomIndexB;

			// range [0 .. keepPopulation - 1]
			randomIndexA = rand() % keepPopulation;
			do {
				randomIndexB = rand() % keepPopulation;
			} while (randomIndexA == randomIndexB);

			population[i] = crossLink(population[randomIndexA], population[randomIndexB]);
		}

		// Mutation time :D
		for(int i = 0; i < numMutations; ++i){

			unsigned randomIndex;

			randomIndex = 1 + rand() % (populationSize - 1);
			population[randomIndex].mutate();
		}

		if (gen % 10 == 0) {
        	cout << "Generation " << gen / 10 << ":  shortest path found "
            	<< population[0].getCircuitLength() << endl;
    	}

	}

	return population[0];
}
/**
 * @description: compare two genome's length, if GA>GB return False, else True.
 * @param {TSPGenome} &GA
 * @param {TSPGenome} &GB
 * @return {*}
 */
bool isShorterPath(const TSPGenome &GA, const TSPGenome &GB) {
	double fitnessGA = GA.getCircuitLength();
	double fitnessGB = GB.getCircuitLength();
	if(fitnessGA < fitnessGB)
		return true;
	else 
		return false; 
}
