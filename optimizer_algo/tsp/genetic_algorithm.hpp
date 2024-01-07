// This Genetic Algorithm class guides the evolution process of DNAs.
// It manages the DNAs of the current population and next generation
// and implements a strategy to select parent DNAs for reproduction.
// It can be used for different problems without modification.


#pragma once

#include <iostream>
using std::cout;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <memory>
using std::unique_ptr;

#include <utility>
using std::move;
using std::pair;
using std::make_pair;

#include <limits>
#include <algorithm>

#include "randnumgen.hpp"


// The objective of an optimization algorithm is to minimize or maximize
// the output value of an objective function. The algorithm has to
// systematically pick or create input variables (solutions) for the
// objective function until a sufficient or optimal solution is found.
//
// When the output value of an objective function should be minimized,
// it is also called a 'loss function' or 'cost function'. If the output
// value should be maximized, the objective function is often called a
// 'reward function' or a 'fitness function' in the case of
// a Genetic Algorithm.

enum GeneticAlgorithmObjective { Minimization, Maximization };


template<typename DNA>
class GeneticAlgorithm
{
public:
    GeneticAlgorithm(
        const DNA* initialDNA,
        GeneticAlgorithmObjective objective,
        size_t populationSize,
        double percentOfBestCanReproduce,
        double recombinationProbability,
        double mutationProbability);

    void evolveNextGeneration();

    GeneticAlgorithmObjective getObjective() const;
    size_t getPopulationSize() const;
    size_t getNumEvolvedGenerations() const;
    size_t getNumFitnessEvaluations() const;
    size_t getNumFitnessImprovements() const;
    size_t getNumGenerationImprovements() const;
    double getBestDNAfitness() const;
    string getBestDNAstring() const;
    bool   isProblemSolved() const;

private:
    inline void createInitialPopulation();
    inline void calcAllFitnessValues();
    inline void createSelectionPool();
    inline void selectParentsFromPool();

private:
    const DNA* _initialDNA;
    GeneticAlgorithmObjective _objective;
    size_t _populationSize;
    size_t _halfPopulationSize;
    double _percentOfBestCanReproduce;
    double _recombinationProbability;
    double _mutationProbability;
    vector<unique_ptr<DNA>> _population;
    vector<unique_ptr<DNA>> _nextGeneration;
    vector<size_t> _selectionPool;
    DNA* _parentDNA1;
    DNA* _parentDNA2;
    DNA* _bestDNA;
    size_t _numEvolvedGenerations;
    size_t _numFitnessEvaluations;
    size_t _numFitnessImprovements;
    size_t _numGenerationImprovements;
};


// ----------------------------------------------------------------------------


template<typename DNA>
GeneticAlgorithm<DNA>::GeneticAlgorithm(
    const DNA* initialDNA,
    GeneticAlgorithmObjective objective,
    size_t populationSize,
    double percentOfBestCanReproduce,
    double recombinationProbability,
    double mutationProbability)
{
    _initialDNA = initialDNA;
    _objective = objective;

    // Ensure a minimum and even population size
    _populationSize = std::max<size_t>(4, populationSize);

    if (_populationSize % 2 != 0)
        _populationSize++;

    _halfPopulationSize = _populationSize / 2;

    if (_populationSize != populationSize)
        cout << "Adjusted population size to an even size: "
             << _populationSize << "\n";

    _percentOfBestCanReproduce = percentOfBestCanReproduce;
    _recombinationProbability = recombinationProbability;
    _mutationProbability = mutationProbability;

    _selectionPool.reserve(_populationSize);

    _parentDNA1 = nullptr;
    _parentDNA2 = nullptr;
    _bestDNA = nullptr;

    _numEvolvedGenerations = 0;
    _numFitnessEvaluations = 0;
    _numFitnessImprovements = 0;
    _numGenerationImprovements = 0;

    createInitialPopulation();
}


template<typename DNA>
GeneticAlgorithmObjective GeneticAlgorithm<DNA>::getObjective() const
{
    return _objective;
}


template<typename DNA>
size_t GeneticAlgorithm<DNA>::getPopulationSize() const
{
    return _populationSize;
}


template<typename DNA>
size_t GeneticAlgorithm<DNA>::getNumEvolvedGenerations() const
{
    return _numEvolvedGenerations;
}


template<typename DNA>
size_t GeneticAlgorithm<DNA>::getNumFitnessEvaluations() const
{
    return _numFitnessEvaluations;
}


template<typename DNA>
size_t GeneticAlgorithm<DNA>::getNumFitnessImprovements() const
{
    return _numFitnessImprovements;
}


template<typename DNA>
size_t GeneticAlgorithm<DNA>::getNumGenerationImprovements() const
{
    return _numGenerationImprovements;
}


template<typename DNA>
double GeneticAlgorithm<DNA>::getBestDNAfitness() const
{
    return _bestDNA->getFitness();
}


template<typename DNA>
string GeneticAlgorithm<DNA>::getBestDNAstring() const
{
    return _bestDNA->toString();
}


template<typename DNA>
bool GeneticAlgorithm<DNA>::isProblemSolved() const
{
    return _bestDNA->isSolved();
}


template<typename DNA>
inline void GeneticAlgorithm<DNA>::createInitialPopulation()
{
    // INITIALIZATION  (See Tutorial)
    _population.resize(_populationSize);
    _nextGeneration.resize(_populationSize);

    for (size_t iDNA = 0; iDNA < _populationSize; iDNA++)
    {
        _population[iDNA] = unique_ptr<DNA>( new DNA(*_initialDNA) );
        _population[iDNA]->initGenesWithRandomValues();

        _nextGeneration[iDNA] = unique_ptr<DNA>( new DNA(*_initialDNA) );
        _nextGeneration[iDNA]->initGenesWithRandomValues();
    }

    // FITNESS EVALUATION  (See Tutorial)
    calcAllFitnessValues();

    _numEvolvedGenerations = 1;
}


template<typename DNA>
inline void GeneticAlgorithm<DNA>::calcAllFitnessValues()
{
    double minFitness = std::numeric_limits<double>::max();
    double maxFitness = std::numeric_limits<double>::min();
    double fitness = 0.0;

    if (_objective == GeneticAlgorithmObjective::Minimization)
    {
        for (size_t iDNA = 0; iDNA < _populationSize; iDNA++)
        {
            fitness = _population[iDNA]->calcFitness();
            _numFitnessEvaluations++;

            if (fitness < minFitness)
            {
                minFitness = fitness;
                _bestDNA = _population[iDNA].get();
                _numFitnessImprovements++;
            }
        }
    }
    else if (_objective == GeneticAlgorithmObjective::Maximization)
    {
        for (size_t iDNA = 0; iDNA < _populationSize; iDNA++)
        {
            fitness = _population[iDNA]->calcFitness();
            _numFitnessEvaluations++;

            if (fitness > maxFitness)
            {
                maxFitness = fitness;
                _bestDNA = _population[iDNA].get();
                _numFitnessImprovements++;
            }
        }
    }
}


template<typename DNA>
void GeneticAlgorithm<DNA>::evolveNextGeneration()
{
    // PRE-SELECTION
    // Potential parent DNAs for reproduction will be picked from the current
    // population and stored in a selection pool. DNAs with better fitness
    // values have a higher chance to be part of the pool. 

    createSelectionPool();

    // In this implementation parents always create two children and the
    // size of the next generation is equal to the current population size.
    // You can experiment with a variable number of children or let the
    // population size increase or decrease over time.

    for (size_t iDNA = 0; iDNA < _halfPopulationSize; iDNA++)
    {
        // SELECTION   (See Tutorial)
        selectParentsFromPool();

        // RECOMBINATION  (See Tutorial)
        auto child1 = _nextGeneration[iDNA].get();
        auto child2 = _nextGeneration[iDNA + _halfPopulationSize].get();

        if (getRandomTrueWithProbability(_recombinationProbability))
        {
            child1->recombineGenes(*_parentDNA1, *_parentDNA2);
            child2->recombineGenes(*_parentDNA2, *_parentDNA1);
        }
        else
        {
            child1->copyGenes(*_parentDNA1);
            child2->copyGenes(*_parentDNA2);
        }
        
        // MUTATION  (See Tutorial)
        child1->mutateGenes(_mutationProbability);
        child2->mutateGenes(_mutationProbability);
    }

    // REPLACEMENT  (See Tutorial)
    _population.swap(_nextGeneration);

    // FITNESS EVALUATION  (See Tutorial)
    const double previousBestFitness = _bestDNA->getFitness();
    calcAllFitnessValues();

    // PERFORMANCE STATISTICS
    if (_bestDNA->getFitness() < previousBestFitness &&
        _objective == GeneticAlgorithmObjective::Minimization)
    {
        _numGenerationImprovements++;
    }
    else if (_bestDNA->getFitness() > previousBestFitness &&
        _objective == GeneticAlgorithmObjective::Maximization)
    {
        _numGenerationImprovements++;
    }

    _numEvolvedGenerations++;
}


template<typename DNA>
inline void GeneticAlgorithm<DNA>::createSelectionPool()
{
    // Usual methods for SELECTION in a Genetic Algorithm are:
    //
    // - Fitness Proportionate Selection
    // - Roulette Wheel Selection
    // - Stochastic Universal Sampling
    // - Rank Selection
    // - Elitist Selection
    // - Tournament Selection
    // - Truncation Selection
    //
    // This implementation combines some of the ideas.
    //
    // When a new generation will be evolved, not every DNA in the population
    // will get a chance to reproduce. A given percentage of the best DNAs will
    // be placed in a pool from which parents will be randomly selected 
    // for reproduction. A DNA with a better fitness value will be inserted
    // more often into the selection pool than a DNA with a low fitness value.
    // A DNA which is represented more often in the pool has a higher
    // probability to be selected as a parent.
    //
    // Imagine the selection pool as a roulette wheel. Each DNA in the pool gets
    // a slot in a roulette wheel. The size of the slot is based on the fitness
    // of the DNA. When selecting parents for reproduction it is like throwing
    // two balls into a roulette wheel. The balls have a higher chance of
    // landing in bigger slots than in smaller slots.
    //
    // When a single DNA has a very high fitness value and all other DNAs have
    // very low values, the one good DNA will become overrepresented in the
    // pool and it will be selected for reproduction nearly every time.
    // This (elitism) will reduce the diversity in the next generation and the
    // evolution may get stuck very fast. To get a more balanced representation
    // of good DNAs in the pool, their relative fitness values (to each other)
    // will be used instead of their absolute fitness values. All DNAs of a
    // population will be sorted in a list based on their fitness values.
    // Their position in the list will define their rank (relative fitness).
    // A percentage of the best ranked DNAs will be inserted into the
    // selection pool.
    //
    // The DNAs will not be copied to the ranking list or the selection pool.
    // Both are vectors of indices pointing to the DNAs in the population.


    // A ranking of DNA indices will be created. The DNAs will be sorted based
    // on their fitness values.

    vector<pair<size_t, double>> ranking(_populationSize);

    for (size_t iDNA = 0; iDNA < _populationSize; iDNA++)
        ranking[iDNA] = make_pair(iDNA, _population[iDNA]->getFitness());

    if (_objective == GeneticAlgorithmObjective::Minimization)
    {
        std::sort(
            ranking.begin(),
            ranking.end(),
            [](pair<size_t, double>& left, pair<size_t, double>& right)
            {
                return left.second < right.second;
            });
    }
    else if (_objective == GeneticAlgorithmObjective::Maximization)
    {
        std::sort(
            ranking.begin(),
            ranking.end(),
            [](pair<size_t, double>& left, pair<size_t, double>& right)
            {
                return left.second > right.second;
            });
    }


    // Take some percentage of the best DNAs from the ranking as candidates
    // for the selection pool. Better ranked candidates have a higher
    // probability to get inserted into the pool.

    _selectionPool.clear();

    const size_t candidateCount =
        static_cast<size_t>( ranking.size() * (_percentOfBestCanReproduce * 0.01) );

    const double candidateQuotient = 1.0 / candidateCount;

    double insertionProbability = 0.0;


    for (size_t rank = 0; rank < candidateCount; rank++)
    {
        insertionProbability = 1.0 - rank * candidateQuotient;

        // Optional: Strong (exponential) falloff
        insertionProbability = std::pow(insertionProbability, 2.0);

        if (getRandomTrueWithProbability(insertionProbability))
            _selectionPool.push_back(ranking[rank].first);
    }
}


template<typename DNA>
inline void GeneticAlgorithm<DNA>::selectParentsFromPool()
{
    // Randomly select two DNAs from the selection pool.
    // DNAs with better fitness values have a higher chance to get selected
    // because they have a higher appearance in the selection pool.
    // It is tolerable that sometimes both parents will be the same DNA.

    const auto poolIdx1 =
        getRandomIntegerInRange<size_t>(0, _selectionPool.size() - 1);

    const auto poolIdx2 =
        getRandomIntegerInRange<size_t>(0, _selectionPool.size() - 1);

    const auto dnaIdx1 = _selectionPool[poolIdx1];
    const auto dnaIdx2 = _selectionPool[poolIdx2];

    _parentDNA1 = _population[dnaIdx1].get();
    _parentDNA2 = _population[dnaIdx2].get();
}