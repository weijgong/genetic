/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-02 01:33:21
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-20 00:38:26
 * @FilePath: /root/genetic/sat_algorithm/main.c
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
// I.E.
#include "sat_algorithm.h"

extern struct TimeWindow tw_list [MAX_TARGET_NUM];

int main(){
    Sense_mode *SenseModeArray = NULL;
    int numMode = 0;
    read_sense_mode(&SenseModeArray,&numMode);
    
    for (int i = 0; i < numMode; i++) {
        printf("Mode %d:\n", i + 1);
        printf("Imaging Duration: %.2lf seconds\n", SenseModeArray[i].imagingDuration);
        printf("Resolution: %.2lf meters\n", SenseModeArray[i].resolution);
        printf("Swath Width: %.2lf kilometers\n", SenseModeArray[i].swathWidth);
        printf("Is Slant View: %d\n", SenseModeArray[i].isSlantView);
        printf("Slant Angle: %.2lf degrees\n", SenseModeArray[i].slantAngle);
        printf("\n");
    }    

    init_time_windows();
    for(int i = 0;i<MAX_TARGET_NUM;i++){
        printf("Target no %d:\n", tw_list[i].target_no);
        printf("Start time windows: %ld \n", tw_list[i].start_time);
        printf("End time windows: %ld \n", tw_list[i].stop_time);
        printf("Time windows duration: %ld \n", tw_list[i].durations);
        printf("\n");
    }

    init_position();
    for(int i = 0;i<MAX_TARGET_NUM;i++){
        printf("Target Latitude: %f \n", pos_list[i].latitude);
        printf("Target Longitude: %f \n", pos_list[i].longitude);
        printf("\n");
    }
    
    free(SenseModeArray);
    return 0;
}