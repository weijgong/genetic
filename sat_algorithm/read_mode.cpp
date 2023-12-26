/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-19 21:52:49
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-26 15:38:41
 * @FilePath: /root/genetic/sat_algorithm/read_mode.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "sat_algorithm.h"

void read_sense_mode(Sense_mode **dataArray,int *numEntries) {
    FILE *file = fopen("/root/genetic/sat_data/sat_mode.csv", "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(1);
    }

    char line[MAX_LINE_LENGTH];
    *numEntries = 0;  

    // 获取模式数
    while (fgets(line, sizeof(line), file) != NULL) {
        (*numEntries)++;
    }

    *dataArray = (Sense_mode *)malloc((*numEntries) * sizeof(Sense_mode));
    if (dataArray == NULL) {
        perror("Memory allocation error");
        fclose(file);
        exit(1);
    }

    // 获取好文件总行数之后将文件指针回到文件开头。
    if (fseek(file, 0, SEEK_SET) != 0) {
        perror("Error resetting file pointer");
        fclose(file);
        exit(1);
    }
    for (int i = 0; i < (*numEntries); i++) {
        fgets(line, sizeof(line), file);
        sscanf(line, "%lf %lf %lf %d %lf",
               &(*dataArray)[i].imagingDuration,
               &(*dataArray)[i].resolution,
               &(*dataArray)[i].swathWidth,
               &(*dataArray)[i].isSlantView,
               &(*dataArray)[i].slantAngle);
    }

    fclose(file);
}