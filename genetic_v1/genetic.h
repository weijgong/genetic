/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2024-01-02 00:21:19
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2024-01-02 12:01:49
 * @FilePath: /gongweijing/genetic/sat_algorithm/genetic.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%A
 */
#ifndef GENETIC_H
#define GENETIC_H
#include "coding.h"
#include <vector>
#include <algorithm>

using namespace std;

struct Individual {
    EvaluationCode* genes;
    double fitness;
};

class GeneticAlgorithm {
private:
    int populationSize;    // 种群大小
    int geneSize;          // 每个个体的基因长度
    double mutationRate;   // 变异率
    double crossoverRate;  // 交叉率
public:
    
    EvaluationCode* getRandomGene(Sense_mode *SenseModeArray);
    GeneticAlgorithm(int popSize, double mutationRate, double crossoverRate);
    vector<Individual> initializePopulation(Sense_mode *SenseModeArray);
    void calculateFitness(Individual& individual);
    void nout_individual();
};
#endif