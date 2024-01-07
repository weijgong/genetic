/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-02 01:33:21
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2024-01-07 10:47:27
 * @FilePath: /gongweijing/genetic/sat_algorithm/main.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
// I.E.
#include "sat_algorithm.h"
#include "coding.h"
#include "genetic.h"

using namespace std;


struct TimeWindow tw_list [MAX_TARGET_NUM];
struct Position   pos_list[MAX_TARGET_NUM];
struct SteoCord   cor_list[MAX_TARGET_NUM];

#include <stdio.h>


int main(){
    srand((unsigned)time(NULL));
    Sense_mode *SenseModeArray = NULL;
    int numMode = 0;
    read_sense_mode(&SenseModeArray,&numMode);
    
    // for (int i = 0; i < numMode; i++) {
    //     printf("Mode %d:\n", i + 1);
    //     printf("Imaging Duration: %.2lf seconds\n", SenseModeArray[i].imagingDuration);
    //     printf("Resolution: %.2lf meters\n", SenseModeArray[i].resolution);
    //     printf("Swath Width: %.2lf kilometers\n", SenseModeArray[i].swathWidth);
    //     printf("Is Slant View: %d\n", SenseModeArray[i].isSlantView);
    //     printf("Slant Angle: %.2lf degrees\n", SenseModeArray[i].slantAngle);
    //     printf("\n");
    // }    

    init_time_windows();
    // for(int i = 0;i<MAX_TARGET_NUM;i++){
    //     printf("Target no %d:\n", tw_list[i].target_no);
    //     printf("Start time windows: %ld \n", tw_list[i].start_time);
    //     printf("End time windows: %ld \n", tw_list[i].stop_time);
    //     printf("Time windows duration: %ld \n", tw_list[i].durations);
    //     printf("\n");
    // }

    extend_time_windows();
    time_t earliest_time_start = (time_t)INFINITY;
    time_t slowes_time_stop = 0;
    for(int i = 0;i<MAX_TARGET_NUM;i++){
        if(earliest_time_start>tw_list[i].start_time){
            earliest_time_start = tw_list[i].start_time;            
        }
        if(slowes_time_stop<tw_list[i].stop_time){
            slowes_time_stop = tw_list[i].stop_time;
        }
        // printf("Target no %d:\n", tw_list[i].target_no);
        // printf("Start time windows: %ld \n", tw_list[i].start_time);
        // printf("End time windows: %ld \n", tw_list[i].stop_time);
        // printf("Time windows duration: %ld \n", tw_list[i].durations);
        // printf("\n");
    }
    // cout<<"最早的时间窗口开始时间为："<<earliest_time_start<<endl;
    // cout<<"最晚的时间窗口结束时间为："<<slowes_time_stop<<endl;

    // 主要作用是进行约束条件和编码后生成开始窗口
    for(int i = 0;i<MAX_TARGET_NUM;i++){
        tw_list[i].start_time -= earliest_time_start;
        tw_list[i].stop_time  -= earliest_time_start;
        // printf("Start:%ld,end:%ld,duration:%ld\n",tw_list[i].start_time,tw_list[i].stop_time,tw_list[i].durations);
    }

    
    // 绘制各个任务同一个timeline的甘特图
    // plot_target_windows("target_time_windows.png",earliest_time_start,slowes_time_stop);

    // init_position();
    // for(int i = 0;i<MAX_TARGET_NUM;i++){
    //     printf("Target Latitude: %f \n", pos_list[i].latitude);
    //     printf("Target Longitude: %f \n", pos_list[i].longitude);
    //     printf("\n");
    // }
    GeneticAlgorithm algo(100,0.1,0.9);
    vector<Individual> population = algo.initializePopulation(SenseModeArray);

    
    // plot_individual(population[0],earliest_time_start,slowes_time_stop);

    algo.assginFitness(population);
    // 输出当前的个体数
    // algo.abort_population(population);

    vector<Individual> newPopulation;
    for(int i = 0;i < ;i ++){
        vector<Individual> parent = algo.TournamentSelection(population,2,6);
        // 输出生成的两个父代,便于进行交叉重组
        // algo.abort_population(parent);

        vector<Individual> child = algo.crossover(parent);
        // cout<<"子代个体交叉后:"<<endl;
        // algo.abort_population(child);

        algo.mutate(child[0]);
        algo.mutate(child[1]);
        // cout<<"子代个体变异后:"<<endl;
        // algo.abort_population(child);

        algo.reGenerateOfferingInfo(child[0]);
        algo.reGenerateOfferingInfo(child[1]);

        algo.assginFitness(child);

        Individual best_child = algo.Tournament(child);
        // algo.abort_population(child);
        // algo.nout_individual(best_child);
        
        newPopulation.push_back(best_child);
    }    
    population = newPopulation;
    // 求解最好的个体就对所有的进行锦标赛，找出其中最好的个体
    Individual best_ind = algo.Tournament(population);
    // cout<<"最好的个体："<<endl;
    // algo.nout_individual(best_ind);
    
    free(SenseModeArray);
    return 0;
}

