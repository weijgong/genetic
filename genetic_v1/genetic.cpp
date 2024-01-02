/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2024-01-02 00:21:15
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2024-01-02 00:41:49
 * @FilePath: /root/genetic/sat_algorithm/genetic.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>


// 定义个体结构
struct Individual {
    std::vector<double> genes;  // 基因
    double fitness;            // 适应度
};

class GeneticAlgorithm {
private:
    int populationSize;    // 种群大小
    int geneSize;          // 每个个体的基因长度
    double mutationRate;   // 变异率
    double crossoverRate;  // 交叉率

public:
    GeneticAlgorithm(int popSize, int genSize, double mutationRate, double crossoverRate)
        : populationSize(popSize), geneSize(genSize), mutationRate(mutationRate), crossoverRate(crossoverRate) {}

    // 初始化种群
    std::vector<Individual> initializePopulation() {
        std::vector<Individual> population;
        for (int i = 0; i < populationSize; ++i) {
            Individual individual;
            individual.genes.resize(geneSize);
            for (int j = 0; j < geneSize; ++j) {
                individual.genes[j] = getRandomGene();
            }
            population.push_back(individual);
        }
        return population;
    }

    // 计算个体适应度
    void calculateFitness(Individual& individual) {
        // 这里使用一个虚构的适应度函数，实际问题需要根据具体情况定义
        individual.fitness = 0.0;
        for (double gene : individual.genes) {
            individual.fitness += gene;
        }
    }

    // 选择
    Individual selectParent(const std::vector<Individual>& population) {
        // 这里使用轮盘赌选择，你可以根据需要选择其他选择算法
        double totalFitness = 0.0;
        for (const Individual& individual : population) {
            totalFitness += individual.fitness;
        }

        double targetFitness = getRandomNumber(0.0, totalFitness);
        double currentFitness = 0.0;

        for (const Individual& individual : population) {
            currentFitness += individual.fitness;
            if (currentFitness >= targetFitness) {
                return individual;
            }
        }

        // 如果没有选择到任何个体，返回最后一个个体
        return population.back();
    }

    // 交叉
    Individual crossover(const Individual& parent1, const Individual& parent2) {
        // 一点交叉，你可以根据需要选择其他交叉算法
        int crossoverPoint = getRandomNumber(1, geneSize - 1);

        Individual child;
        child.genes.resize(geneSize);

        for (int i = 0; i < crossoverPoint; ++i) {
            child.genes[i] = parent1.genes[i];
        }

        for (int i = crossoverPoint; i < geneSize; ++i) {
            child.genes[i] = parent2.genes[i];
        }

        return child;
    }

    // 变异
    void mutate(Individual& individual) {
        for (double& gene : individual.genes) {
            if (getRandomNumber(0.0, 1.0) < mutationRate) {
                gene = getRandomGene();
            }
        }
    }

    // 遗传算法主循环
    Individual evolve(int generations) {
        std::vector<Individual> population = initializePopulation();

        for (int generation = 0; generation < generations; ++generation) {
            for (Individual& individual : population) {
                calculateFitness(individual);
            }

            std::vector<Individual> newPopulation;

            for (int i = 0; i < populationSize; ++i) {
                Individual parent1 = selectParent(population);
                Individual parent2 = selectParent(population);

                Individual child = crossover(parent1, parent2);
                mutate(child);

                newPopulation.push_back(child);
            }

            population = newPopulation;
        }

        // 返回最终种群中适应度最高的个体
        return *std::max_element(population.begin(), population.end(),
                                  [](const Individual& ind1, const Individual& ind2) {
                                      return ind1.fitness < ind2.fitness;
                                  });
    }

    // 获取随机基因
    double getRandomGene() {
        return getRandomNumber(-1.0, 1.0);
    }

    // 获取指定范围内的随机数
    double getRandomNumber(double min, double max) {
        return min + static_cast<double>(rand()) / RAND_MAX * (max - min);
    }
};

int main() {
    // 设置随机种子
    srand(static_cast<unsigned>(time(nullptr)));

    // 创建遗传算法对象
    GeneticAlgorithm geneticAlgorithm(50, 10, 0.01, 0.8);

    // 运行遗传算法，进行100代演化
    Individual result = geneticAlgorithm.evolve(100);

    // 打印最终结果
    std::cout << "最优个体的适应度: " << result.fitness << std::endl;
    std::cout << "最优个体的基因: ";
    for (double gene : result.genes) {
        std::cout << gene << " ";
    }
    std::cout << std::endl;

    return 0;
}