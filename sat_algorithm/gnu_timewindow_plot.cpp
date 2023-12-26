/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-26 15:34:19
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-26 15:34:32
 * @FilePath: /root/genetic/sat_algorithm/gnu_timewindow_plot.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "sat_algorithm.h"

void generateDataFile(const char* filename) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file\n");
        return;
    }

    for (int i = 0; i < MAX_TARGET_NUM; ++i) {
        fprintf(file, "%d %d\n", tw_list[i].start_time, tw_list[i].stop_time);
    }

    fclose(file);
}