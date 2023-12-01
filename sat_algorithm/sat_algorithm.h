/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-02 01:33:32
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-02 01:49:09
 * @FilePath: /root/genetic/sat_algorithm/sat_algorithm.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <stdio.h>
#include <time.h>

#define MAX_LINE_LENGTH 256


// 定义结构体表示闰秒的时间点和秒数
struct LeapSecond {
    time_t time;
    int seconds;
};

struct TimeWindow{
    time_t start_time;
    time_t stop_time;
    time_t durations;
};

// 注意转换成世界秒之后再进行间隔的计算即可。
struct AccessRecord {
    // 序号
    int access;
    // 开始时间
    char start_time[MAX_LINE_LENGTH];
    // 结束时间
    char stop_time[MAX_LINE_LENGTH];
    // 间隔时间
    double duration;
};

time_t utc_to_tai(const char *utc_time);
