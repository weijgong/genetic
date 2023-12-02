/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-02 01:33:15
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-03 00:41:46
 * @FilePath: /root/genetic/sat_algorithm/read_access.c
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "sat_algorithm.h"


void parse_csv_line(char *line, struct AccessRecord *record) {
    char *token;
    token = strtok(line, ",");
    record->target_no = atoi(token);

    token = strtok(NULL, ",");
    strncpy(record->start_time, token, MAX_LINE_LENGTH);

    token = strtok(NULL, ",");
    strncpy(record->stop_time, token, MAX_LINE_LENGTH);
    
    token = strtok(NULL,",");
    record->duration = atoi(token);
}

void init_time_windows(){
        // 文件指针
    FILE *file;

    // 文件路径
    const char *file_path = "/root/genetic/sat_data/Satellite-Satellite1-Sensor-Sensor1-To-Target-Target1_Access.csv";

    // 打开文件以供读取
    file = fopen(file_path, "r");

    // 检查文件是否成功打开
    if (file == NULL) {
        fprintf(stderr, "无法打开文件：%s\n", file_path);
        return 1;  // 返回非零值表示程序出错
    }

    char line[MAX_LINE_LENGTH];

    int n = 0;
    // 读取CSV文件
    while (fgets(line, sizeof(line), file) != NULL) {
        struct AccessRecord record;
        parse_csv_line(line, &record);
        // printf("%d %s %s\n",record.target_no,record.start_time,record.stop_time);
        tw_list[n].start_time=utc_to_tai(record.start_time);
        tw_list[n].stop_time=utc_to_tai(record.stop_time);
        tw_list[n].durations=tw_list[n].stop_time-tw_list[n].start_time;
        // printf("%ld %ld\n",utc_to_tai(record.start_time),utc_to_tai(record.stop_time));
        n++;
    }
    // 关闭文件
    fclose(file);
}

// int main() {
//     init_time_windows();
//     return 0;  // 返回零表示程序执行成功
// }
