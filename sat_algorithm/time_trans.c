/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-01 23:39:37
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-02 00:54:07
 * @FilePath: /root/genetic/sat_algorithm/time_trans.c
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <stdio.h>
#include <time.h>

// 定义结构体表示闰秒的时间点和秒数
struct LeapSecond {
    time_t time;
    int seconds;
};

// 函数：将UTC转换为TAI
time_t utc_to_tai(const char *utc_time) {
    // 已知的闰秒表
    struct LeapSecond leap_seconds[] = {
        {1288348800, 16},  // 2012-07-01 00:00:00
        // 添加其他已知的闰秒时间
    };

    int num_leap_seconds = sizeof(leap_seconds) / sizeof(leap_seconds[0]);

    struct tm tm_struct = {0};
    time_t utc_seconds;

    // 解析UTC时间字符串
    if (strptime(utc_time, "%Y-%m-%d %H:%M:%S", &tm_struct) == NULL) {
        fprintf(stderr, "解析日期时间失败\n");
        return -1;
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
// I.E.
int main() {
    // 示例UTC时间
    const char *utc_time = "2023-01-01 12:00:00";

    // 转换为TAI
    time_t tai_time = utc_to_tai(utc_time);

    // 打印结果
    printf("UTC时间: %s\n", utc_time);
    printf("TAI时间: %s", asctime(gmtime(&tai_time)));  // 使用gmtime将time_t转换为struct tm
    printf("second is:%ld\n",tai_time);
    return 0;
}

