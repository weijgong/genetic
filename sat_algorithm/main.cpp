/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-02 01:33:21
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2024-01-02 00:38:18
 * @FilePath: /root/genetic/sat_algorithm/main.c
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
// I.E.
#include "sat_algorithm.h"
#include "coding.h"

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

    EvaluationCode* ec_simple = (EvaluationCode*)malloc(sizeof(EvaluationCode)*MAX_TARGET_NUM);
    for(int i = 0;i < MAX_TARGET_NUM; i++){
        // 种群初始化
        ec_simple[i].init_individual(SenseModeArray,i);
    }


    free(SenseModeArray);
    return 0;
}

