/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-02 01:33:32
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2024-01-08 17:05:54
 * @FilePath: /gongweijing/genetic/sat_algorithm/sat_algorithm.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef SAT_ALGORITHM_H

#define SAT_ALGORITHM_H

#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <bits/stdc++.h>
#include <unistd.h>

using namespace std;

#define MAX_LINE_LENGTH 256
#define MAX_TARGET_NUM 5

# define M_PI                3.14159265358979323846	
# define EARTH_RADIUS_METER  6371000
# define INCLINE_ANGLE       -35*M_PI/180

// 定义结构体表示闰秒的时间点和秒数
struct LeapSecond {
    time_t time;
    int seconds;
};

struct TimeWindow{
    int target_no;
    time_t start_time;
    time_t stop_time;
    time_t durations;
};

// 经纬度坐标
struct Position{
    // 纬度
    float latitude;
    // 经度
    float longitude;
};

// 转换为三维坐标系
struct SteoCord{
    float x;
    float y;
    float z;
};

struct Distance{
    float along_orbit;
    float vertical_orbit;
    float true_dist;
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


typedef struct Sense_mode{
    double imagingDuration;//成像时长(秒)
    double resolution;     //分辨率(m)
    double swathWidth;     //幅宽(km)
    int isSlantView;       //是否可以斜视(1是/0否)
    double slantAngle;     //斜视角度(deg)
} Sense_mode;

extern struct TimeWindow tw_list [MAX_TARGET_NUM];
extern struct Position   pos_list[MAX_TARGET_NUM];
extern struct SteoCord   cor_list[MAX_TARGET_NUM];

typedef struct Distance Distance;
typedef struct TimeWindow TimeWindow;
typedef struct Position Position;
typedef struct SteoCord SteoCord;

// 获取工作目录
string getExecuatePath();
// 获取工作目录上级目录
string getExecuateParentPath();

// 将UTC转换为TAI
time_t utc_to_tai(char *utc_time);
// 读取时间窗口，并将其赋值到tw_list中
void init_time_windows();
// 读取时间窗口的文件
void parse_csv_line(char *line, struct AccessRecord *record);
// 将时间窗口进行长度扩展，便于后续进行调度
void extend_time_windows();
// 写入时间窗口dat文件
void generateDataFile(const char* filename);
// 绘制时间窗口
void plot_time_window(const char* filename);
// 绘制各个目标的时间窗口
void plot_target_windows(const char* filename,int earliest_time_start,int slowes_time_stop);

// 读取经纬度文件，并将其赋值到tw_list中
void init_position();
// 角度转为弧度
double to_radians(double degrees);
// 弧度转为角度
double to_degree(double radians);
// 根据haversine求解的球面距离
double compute_haversine_distance(struct Position x, struct Position y);
// 根据另一个公式求解球面距离
double compute_curve_distance(struct Position x, struct Position y);
// 计算两点之间的角度
double compute_sphere_angle(struct Position x, struct Position y);
// 经纬度转换为三维坐标系坐标,x轴为0°经度在赤道的投影，z轴为北极方向。
SteoCord to_steocord(struct Position x);
// 三维坐标系的欧氏距离
double compute_steo_dist(SteoCord a,SteoCord b);
// 计算两向量的内积
double compute_steo_product(SteoCord a,SteoCord b);
// 计算两向量差的向量
SteoCord compute_steo_delta_vector(SteoCord b,SteoCord a);
// 计算向量模长
double compute_Vector_modulus(SteoCord vec);
// 计算向量与赤道的线面角
double compute_steo_angle_of_equator(SteoCord vec);
// 计算两向量夹角(deg)
double compute_steo_vectors_angle(SteoCord a,SteoCord b);
// 根据几何定理求夹角
double compute_steo_vector_with_inclination(Position a,Position b);
// 确定最终的角度，以水平方向东北方向为正方向
double compute_steo_real_angle(Position a,Position b);
// 计算投影的长度
Distance compute_project_distance(Position a,Position b);
// 读取传感器观测模式
void read_sense_mode(Sense_mode **dataArray,int *numEntries);


#endif