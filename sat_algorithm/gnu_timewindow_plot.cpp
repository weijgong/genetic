/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-26 15:34:19
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2024-01-07 12:02:07
 * @FilePath: /root/genetic/sat_algorithm/gnu_timewindow_plot.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "sat_algorithm.h"
#include "genetic.h"

void generateDataFile(const char* filename) {
    // FILE* file = fopen(filename, "w");
    // if (file == NULL) {
    //     fprintf(stderr, "Error opening file\n");
    //     return;
    // }

    // for (int i = 0; i < MAX_TARGET_NUM; ++i) {
    //     fprintf(file, "%d %d\n", tw_list[i].start_time, tw_list[i].stop_time);
    // }

    // fclose(file);
}

void plot_time_window(const char* filename){
    // // 调用 gnuplot 命令
    // FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    // if (gnuplotPipe == NULL) {
    //     fprintf(stderr, "Error opening gnuplot\n");
    //     exit(1);
    // }

    // fprintf(gnuplotPipe, "plot '%s' with boxes\n",filename);
    // fprintf(gnuplotPipe, "pause mouse close\n");

    // fclose(gnuplotPipe);
}

void plot_target_windows(const char* filename,int earliest_time_start,int slowes_time_stop){
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

    if (gnuplotPipe != NULL) {
        const char *colors[10] = {
            "red", "green", "blue", "orange", "purple",
            "cyan", "magenta", "yellow", "brown", "pink"
        };

        fprintf(gnuplotPipe, "set terminal pngcairo enhanced font 'arial,10'\n");
        fprintf(gnuplotPipe, "set output '%s'\n",filename);
        // 设定坐标系上下限
        fprintf(gnuplotPipe,"set xrange[%d:%d]\nset yrange[%d:%d]\n",0-1,slowes_time_stop-earliest_time_start,
                0-1,MAX_TARGET_NUM
                );        
        // 设定坐标轴名称
        fprintf(gnuplotPipe,"set xlabel 'timeline'\n set ylabel 'target number'\n");
        // 遍历生成多个矩形
        for(int i = 0;i<MAX_TARGET_NUM;i++){
            fprintf(gnuplotPipe, "set object %d rect from %d,%f to %d,%f fc rgb '%s'\n", 
                    i+1,
                    tw_list[i].start_time,i-0.25,tw_list[i].stop_time,i+0.25,
                    colors[i]
                    );
            // 给每个任务添加注释
            fprintf(gnuplotPipe, "set label '%d' at %f,%d center font ',8'\n", 
                    i+1, (tw_list[i].start_time + tw_list[i].stop_time) / 2.0, i);
        }        
        
        fprintf(gnuplotPipe, "plot '-'\n");
        fprintf(gnuplotPipe, "0 0\n");
        fprintf(gnuplotPipe, "e\n");
        fprintf(gnuplotPipe, "pause mouse keypress 'Press any key to exit...'\n");
        
        fclose(gnuplotPipe);
    } else {
        printf("Error: Couldn't open Gnuplot pipe.\n");
    }
}

void plot_individual(Individual ec,int earliest_time_start,int slowes_time_stop){
    const char* filename = "population_example.png";
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

    if (gnuplotPipe != NULL) {
        const char *colors[10] = {
            "red", "green", "blue", "orange", "purple",
            "cyan", "magenta", "yellow", "brown", "pink"
        };

        fprintf(gnuplotPipe, "set terminal pngcairo enhanced font 'arial,10'\n");
        fprintf(gnuplotPipe, "set output '%s'\n",filename);
        // 设定坐标系上下限
        fprintf(gnuplotPipe,"set xrange[%d:%d]\nset yrange[%d:%d]\n",0-1,slowes_time_stop-earliest_time_start,
                0-1,MAX_TARGET_NUM
                );        
        // 设定坐标轴名称
        fprintf(gnuplotPipe,"set xlabel 'timeline'\n set ylabel 'target number'\n");
        // 遍历生成多个矩形
        for(int i = 0;i<MAX_TARGET_NUM;i++){
            fprintf(gnuplotPipe, "set object %d rect from %f,%f to %f,%f fc rgb '%s'\n", 
                    i+1,
                    ec.genes[i].real_start_time,i-0.25,
                    ec.genes[i].real_finish_time,i+0.25,
                    colors[i]
                    );
            // 给每个任务添加注释
            fprintf(gnuplotPipe, "set label '%d' at %f,%d center font ',8'\n", 
                    i+1, 
                    (ec.genes[i].real_start_time + ec.genes[i].real_finish_time) / 2.0,
                    i);
        }        
        
        fprintf(gnuplotPipe, "plot '-'\n");
        fprintf(gnuplotPipe, "0 0\n");
        fprintf(gnuplotPipe, "e\n");
        fprintf(gnuplotPipe, "pause mouse keypress 'Press any key to exit...'\n");
        
        fclose(gnuplotPipe);
    } else {
        printf("Error: Couldn't open Gnuplot pipe.\n");
    }
}

void plot_individual_with_name(char* filename,Individual ec,int earliest_time_start,int slowes_time_stop){
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");

    if (gnuplotPipe != NULL) {
        const char *colors[10] = {
            "red", "green", "blue", "orange", "purple",
            "cyan", "magenta", "yellow", "brown", "pink"
        };

        fprintf(gnuplotPipe, "set terminal pngcairo enhanced font 'arial,10'\n");
        fprintf(gnuplotPipe, "set output '%s'\n",filename);
        // 设定坐标系上下限
        fprintf(gnuplotPipe,"set xrange[%d:%d]\nset yrange[%d:%d]\n",0-1,slowes_time_stop-earliest_time_start,
                0-1,MAX_TARGET_NUM
                );        
        // 设定坐标轴名称
        fprintf(gnuplotPipe,"set xlabel 'timeline'\n set ylabel 'target number'\n");
        // 遍历生成多个矩形
        for(int i = 0;i<MAX_TARGET_NUM;i++){
            fprintf(gnuplotPipe, "set object %d rect from %f,%f to %f,%f fc rgb '%s'\n", 
                    i+1,
                    ec.genes[i].real_start_time,i-0.25,
                    ec.genes[i].real_finish_time,i+0.25,
                    colors[i]
                    );
            // 给每个任务添加注释
            fprintf(gnuplotPipe, "set label '%d' at %f,%d center font ',8'\n", 
                    i+1, 
                    (ec.genes[i].real_start_time + ec.genes[i].real_finish_time) / 2.0,
                    i);
        }        
        
        fprintf(gnuplotPipe, "plot '-'\n");
        fprintf(gnuplotPipe, "0 0\n");
        fprintf(gnuplotPipe, "e\n");
        fprintf(gnuplotPipe, "pause mouse keypress 'Press any key to exit...'\n");
        
        fclose(gnuplotPipe);
    } else {
        printf("Error: Couldn't open Gnuplot pipe.\n");
    }
}