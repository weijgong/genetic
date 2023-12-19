#include "sat_algorithm.h"

void read_sense_mode(Sense_mode **dataArray,int *numEntries) {
    FILE *file = fopen("/root/genetic/sat_data/sat_mode.csv", "r");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
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
        return 1;
    }

    // 获取好文件总行数之后将文件指针回到文件开头。
    if (fseek(file, 0, SEEK_SET) != 0) {
        perror("Error resetting file pointer");
        fclose(file);
        return 1;
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