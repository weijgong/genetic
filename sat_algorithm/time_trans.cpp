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

time_t tm_to_seconds(string date_) {
    struct tm* timeinfo;
    if (strptime(date_.data(), "%d %b %Y %H:%M:%S", timeinfo) != NULL) {        
    } else {
        printf("Failed to parse date and time.\n");
    }
    return mktime(timeinfo);
}

// 将 struct tm 转换为 Modified Julian Date (MJD)
double tm_to_mjd(string date_) {
    struct tm timeinfo;
    if (strptime(date_.data(), "%d %b %Y %H:%M:%S", &timeinfo) != NULL) {        
    } else {
        printf("Failed to parse date and time.\n");
    }

    int year = timeinfo.tm_year + 1900;
    int month = timeinfo.tm_mon + 1;
    int day = timeinfo.tm_mday;
    int hour = timeinfo.tm_hour;
    int minute = timeinfo.tm_min;
    int second = timeinfo.tm_sec;

    // 如果月份是1月或2月，则年份和月份都要减1
    if (month <= 2) {
        year--;
        month += 12;
    }

    // 计算世纪数和年份的日数
    int a = year / 100;
    int b = 2 - a + a / 4;
    int c = 365.25 * (year + 4716);
    int d = 30.6001 * (month + 1);

    // 计算总天数
    double jd = b + c + d + day - 1524.5;

    // 转换为 Modified Julian Date (MJD)
    double mjd = jd - 2400000.5;

    // 加上小时、分钟和秒的分数部分
    double frac_day = (hour + minute / 60.0 + second / 3600.0) / 24.0;
    mjd += frac_day;

    return mjd;
}

