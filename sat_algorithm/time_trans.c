/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-01 23:39:37
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-02 08:48:18
 * @FilePath: /root/genetic/sat_algorithm/time_trans.c
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <stdio.h>
#include <time.h>
#include "sat_algorithm.h"

// 函数：将UTC转换为TAI
time_t utc_to_tai(char *utc_time) {
    // 已知的闰秒表
    struct LeapSecond leap_seconds[] = {
        {1288348800, 16},  // 2012-07-01 00:00:00
        // 添加其他已知的闰秒时间
    };

    int num_leap_seconds = sizeof(leap_seconds) / sizeof(leap_seconds[0]);

    struct tm tm_struct = {0};
    time_t utc_seconds;
    double mill;

    // 解析UTC时间字符串
    if (sscanf(utc_time, "%d %3s %d %2s:%2s:%2s.%lf", &tm_struct.tm_mday, tm_struct.tm_mon, &tm_struct.tm_year,
                                               &tm_struct.tm_hour, &tm_struct.tm_min, &tm_struct.tm_sec,
                                               &mill) == 7) {
        printf("Day: %d\n", tm_struct.tm_mday);
        printf("Month: %s\n", tm_struct.tm_mon);
        printf("Year: %d\n", tm_struct.tm_year);
    } else {
        printf("解析失败\n");
    }

    // 将解析的时间转换为秒数
    utc_seconds = mktime(&tm_struct);

    // 根据已知的闰秒表，调整时间
    for (int i = 0; i < num_leap_seconds; ++i) {
        if (utc_seconds >= leap_seconds[i].time) {
            utc_seconds += leap_seconds[i].seconds;
        }
    }

    return utc_seconds;
}
// // I.E.
int main() {
    // 示例UTC时间
    char *utc_time = "1 Dec 2023 04:44:17.748";

    // 转换为TAI
    time_t tai_time = utc_to_tai(utc_time);

    // 打印结果
    printf("UTC时间: %s\n", utc_time);
    printf("TAI时间: %s", asctime(gmtime(&tai_time)));  // 使用gmtime将time_t转换为struct tm
    printf("second is:%ld\n",tai_time);
    return 0;
}
