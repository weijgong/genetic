// #include <stdio.h>
// #include <stdlib.h>
// #include <time.h>

// #define POPULATION_SIZE 50
// #define MAIN_CODING_SIZE 5
// #define SUB_CODING_SIZE 5
// #define GENERATIONS 100
// #define CROSSOVER_RATE 0.8
// #define MUTATION_RATE 0.1
// #define TIME_WINDOW_CENTER 10 // 特定时间窗口的中心位置

// // 个体结构体
// typedef struct {
//     int *mainCoding;
//     int *subCoding;
//     int fitness;
// } Individual;

// // 计算适应度函数
// int calculateFitness(Individual *individual) {
//     int combinedValue = 0;
//     for (int i = 0; i < MAIN_CODING_SIZE; ++i) {
//         combinedValue += individual->mainCoding[i] * individual->subCoding[i];
//     }
//     return combinedValue;
// }

// // 生成受约束的子编码
// void generateConstrainedSubCoding(int mainCoding, int *subCoding) {
//     if (mainCoding == 0 || mainCoding == 1) {
//         int center = TIME_WINDOW_CENTER;
//         for (int i = 0; i < SUB_CODING_SIZE; ++i) {
//             subCoding[i] = center + rand() % 7 - 3;
//         }
//     } else if (mainCoding == 2) {
//         for (int i = 0; i < SUB_CODING_SIZE; ++i) {
//             subCoding[i] = rand() % 21;
//         }
//     }
// }

// // 初始化个体
// void initializeIndividual(Individual *individual) {
//     individual->mainCoding = (int *)malloc(MAIN_CODING_SIZE * sizeof(int));
//     individual->subCoding = (int *)malloc(SUB_CODING_SIZE * sizeof(int));

//     for (int j = 0; j < MAIN_CODING_SIZE; ++j) {
//         individual->mainCoding[j] = rand() % 3; // 主编码为0、1、2
//         generateConstrainedSubCoding(individual->mainCoding[j], individual->subCoding);
//     }
//     individual->fitness = calculateFitness(individual);
// }

// // 初始化种群
// void initializePopulation(Individual *population) {
//     for (int i = 0; i < POPULATION_SIZE; ++i) {
//         initializeIndividual(&population[i]);
//     }
// }

// // 选择操作，采用轮盘赌选择
// void selectPopulation(Individual *population, Individual *selectedPopulation) {
//     int totalFitness = 0;
//     for (int i = 0; i < POPULATION_SIZE; ++i) {
//         totalFitness += population[i].fitness;
//     }

//     for (int i = 0; i < POPULATION_SIZE; ++i) {
//         int selected = 0;
//         int targetFitness = rand() % totalFitness;

//         for (int j = 0; j < POPULATION_SIZE; ++j) {
//             selected += population[j].fitness;
//             if (selected >= targetFitness) {
//                 selectedPopulation[i] = population[j];
//                 break;
//             }
//         }
//     }
// }

// // 交叉操作，采用单点交叉
// void crossover(Individual *parent1, Individual *parent2, Individual *child1, Individual *child2) {
//     int crossoverPoint = rand() % MAIN_CODING_SIZE;
//     for (int i = 0; i < crossoverPoint; ++i) {
//         child1->mainCoding[i] = parent1->mainCoding[i];
//         child2->mainCoding[i] = parent2->mainCoding[i];
//         generateConstrainedSubCoding(child1->mainCoding[i], child1->subCoding);
//         generateConstrainedSubCoding(child2->mainCoding[i], child2->subCoding);
//     }
//     for (int i = crossoverPoint; i < MAIN_CODING_SIZE; ++i) {
//         child1->mainCoding[i] = parent2->mainCoding[i];
//         child2->mainCoding[i] = parent1->mainCoding[i];
//         generateConstrainedSubCoding(child1->mainCoding[i], child1->subCoding);
//         generateConstrainedSubCoding(child2->mainCoding[i], child2->subCoding);
//     }

//     child1->fitness = calculateFitness(child1);
//     child2->fitness = calculateFitness(child2);
// }

// // 变异操作，采用单点变异
// void mutate(Individual *individual) {
//     int mutationPoint = rand() % MAIN_CODING_SIZE;
//     individual->mainCoding[mutationPoint] = rand() % 3; // 主编码为0、1、2
//     generateConstrainedSubCoding(individual->mainCoding[mutationPoint], individual->subCoding);
//     individual->fitness = calculateFitness(individual);
// }

// // 释放个体内存
// void freeIndividual(Individual *individual) {
//     free(individual->mainCoding);
//     free(individual->subCoding);
// }

// // 释放种群内存
// void freePopulation(Individual *population) {
//     for (int i = 0; i < POPULATION_SIZE; ++i) {
//         freeIndividual(&population[i]);
//     }
// }

// // 打印个体信息
// void printIndividual(Individual *individual) {
//     printf("Individual: ");
//     for (int j = 0; j < MAIN_CODING_SIZE; ++j) {
//         printf("%d (", individual->mainCoding[j]);
//         for (int k = 0; k < SUB_CODING_SIZE; ++k) {
//             printf("%d ", individual->subCoding[k]);
//         }
//         printf(") ");
//     }
//     printf(" - Fitness: %d\n", individual->fitness);
// }

// // // 主要遗传算法循环
// // int main() {
// //     srand(time(NULL));

// //     Individual *population = (Individual *)malloc(POPULATION_SIZE * sizeof(Individual));
// //     Individual *selectedPopulation = (Individual *)malloc(POPULATION_SIZE * sizeof(Individual));

// //     initializePopulation(population);

// //     for (int generation = 0; generation < GENERATIONS; ++generation) {
// //         selectPopulation(population, selectedPopulation);

// //         for (int i = 0; i < POPULATION_SIZE; i += 2) {
// //             Individual parent1 = selectedPopulation[i];
// //             Individual parent2 = selectedPopulation[i + 1];

// //             Individual child1, child2;
// //             if ((rand() / (double)RAND_MAX) < CROSSOVER_RATE) {
// //                 initializeIndividual(&child1);
// //                 initializeIndividual(&child2);
// //                 crossover(&parent1, &parent2, &child1, &child2);
// //             } else {
// //                 initializeIndividual(&child1);
// //                 initializeIndividual(&child2);
// //             }

// //             if ((rand() / (double)RAND_MAX) < MUTATION_RATE) {
// //                 mutate(&child1);
// //             }
// //             if ((rand() / (double)RAND_MAX) < MUTATION_RATE) {
// //                 mutate(&child2);
// //             }

// //             // 更新种群
// //             freeIndividual(&population[i]);
// //             freeIndividual(&population[i + 1]);
// //             population[i] = child1;
// //             population[i + 1] = child2;
// //         }
// //     }

// //     // 打印最终种群和最优个体
// //     printf("Final Population:\n");
// //     for (int i = 0; i < POPULATION_SIZE; ++i) {
// //         printIndividual(&population[i]);
// //     }

// //     int bestIndex = 0;
// //     for (int i = 1; i < POPULATION_SIZE; ++i) {
// //         if (population[i].fitness > population[bestIndex].fitness) {
// //             bestIndex = i;
// //         }
// //     }

// //     printf("\nFinal Best Individual:\n");
// //     printIndividual(&population[bestIndex]);

// //     // 释放内存
// //     freePopulation(population);
// //     free(selectedPopulation);

// //     return 0;
// // }
