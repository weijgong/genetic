/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-01 23:39:37
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2024-01-06 23:40:04
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

    time_t utc_seconds;

    struct tm tm_time;

    // 使用strptime解析时间字符串
    if (strptime(utc_time, "%d %b %Y %T", &tm_time) == NULL) {
        fprintf(stderr, "解析时间失败\n");
        return 1;
    }

    // 将解析的时间转换为秒数
    utc_seconds = mktime(&tm_time);

    // 根据已知的闰秒表，调整时间
    for (int i = 0; i < num_leap_seconds; ++i) {
        if (utc_seconds >= leap_seconds[i].time) {
            utc_seconds += leap_seconds[i].seconds;
        }
    }

    return utc_seconds;
}

// int main() {
//     char *utc_time = "1 Dec 2023 04:44:17.748";
//     time_t tai = utc_to_tai(utc_time);
//     printf("%ld\n",tai);
//     return 0;
// }

