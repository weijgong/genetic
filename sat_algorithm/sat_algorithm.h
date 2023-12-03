/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-02 01:33:32
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-03 12:08:55
 * @FilePath: /root/genetic/sat_algorithm/sat_algorithm.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <stdio.h>
#include <time.h>

#define MAX_LINE_LENGTH 256
#define MAX_TARGET_NUM 5


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

struct Position{
    // 纬度
    float latitude;
    // 经度
    float longitude;
};

// 注意转换成世界秒之后再进行间隔的计算即可。
struct AccessRecord {
    // 序号
    int target_no;
    // 开始时间
    char start_time[MAX_LINE_LENGTH];
    // 结束时间
    char stop_time[MAX_LINE_LENGTH];
    // 间隔时间
    double duration;
};

struct TimeWindow tw_list[MAX_TARGET_NUM];
struct Position pos_list[MAX_TARGET_NUM];

// 将UTC转换为TAI
time_t utc_to_tai(char *utc_time);
// 读取时间窗口，并将其赋值到tw_list中
void init_time_windows();
// 读取时间窗口的文件
void parse_csv_line(char *line, struct AccessRecord *record);
// 读取经纬度文件，并将其赋值到tw_list中
void init_position();