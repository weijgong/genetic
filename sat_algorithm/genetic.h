/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2024-01-02 00:21:19
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2024-01-07 04:05:32
 * @FilePath: /gongweijing/genetic/sat_algorithm/genetic.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%A
 */
#ifndef GENETIC_H
#define GENETIC_H

#include "sat_algorithm.h"
#include "coding.h"

using namespace std;

typedef struct Individual {
    EvaluationCode* genes;
    double fitness;
}Individual;


// 绘制任意个体的时序图
void plot_individual(Individual ec,int earliest_time_start,int slowes_time_stop);

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
    
    // fitness计算
    double calculateFitness(Individual& individual);
    void assginFitness(vector<Individual>& Population);
    
    // 输出个体跟种群
    void nout_individual(Individual ec);
    void abort_population(vector<Individual>Population);

    // 将错误的个体更正
    Individual RepairUnfeasibleSolution(Individual ec);
    void reGenerateOfferingInfo(Individual& ind);

    // 选择算法
    vector<Individual> SampleTournament(vector<Individual> population,int pool_number);
    Individual Tournament(vector<Individual> TournamentCandidate);
    vector<Individual> TournamentSelection(vector<Individual> population,
                                                         int select_number,
                                                         int tournament_pool_size
                                                         );
    // 交叉算法
    vector<Individual> crossover(vector<Individual> parent);

    // 变异算法
    void mutate(Individual& individual);

    // 主迭代
    
};


#endif