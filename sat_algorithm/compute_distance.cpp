/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-03 21:16:52
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-04 19:10:21
 * @FilePath: /gongweijing/nsga2/sat_algorithm/compute_distance.c
 * @Description: 
 * 
 * Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
 */
#include "sat_algorithm.h"
#include <math.h>

double to_radians(double degrees) {
    return degrees * M_PI / 180.0;
}

double to_degree(double radians) {
    return radians / (M_PI / 180.0);
}

double compute_haversine_distance(struct Position x, struct Position y){
    double d_lat = to_radians(y.latitude - x.latitude);
    double d_lon = to_radians(y.longitude - x.longitude);
    double a = sin(d_lat/2)*sin(d_lat/2) + 
               cos(to_radians(x.latitude))*cos(to_radians(y.latitude))*sin(d_lon/2)*sin(d_lon/2);
    double rad = 2*atan2(sqrt(a),sqrt(1-a));
    double dist = rad*EARTH_RADIUS_METER;
    return dist;
} 

double compute_curve_distance(struct Position x, struct Position y){
    double d_lon = to_radians(y.longitude - x.longitude);
    double a = cos(to_radians(x.latitude))*cos(to_radians(y.latitude))*cos(d_lon)+
               sin(to_radians(x.latitude))*sin(to_radians(y.latitude));
    double rad = acos(a);
    double dist = rad*EARTH_RADIUS_METER;
    return dist;
}

// 计算三维坐标系下面的欧式距离
double compute_steo_dist(SteoCord a,SteoCord b){
    double dist = 0;
    dist = sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
    return dist;
}

// 计算两向量的内积
double compute_steo_product(SteoCord a,SteoCord b){
    double res = a.x*b.x+a.y*b.y+a.z*b.z;
    return res;
}

// 计算两向量差的向量
SteoCord compute_steo_delta_vector(SteoCord b,SteoCord a){
    SteoCord vec;
    vec.x = a.x-b.x;
    vec.y = a.y-b.y;
    vec.z = a.z-b.z;
    return vec;
}

// 计算向量模长
double compute_Vector_modulus(SteoCord vec){
    return sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
}

// 计算向量与赤道的线面角
double compute_steo_angle_of_equator(SteoCord vec){
    SteoCord eq_h={0,0,1};
    double prod = compute_steo_product(eq_h,vec);
    double rad = acos(prod/compute_Vector_modulus(vec));
    return rad;
}

// 计算两向量夹角(rad)
double compute_steo_vectors_angle(SteoCord a,SteoCord b){
    double prod = compute_steo_product(a,b);
    double angle = acos(prod/(compute_Vector_modulus(a)*compute_Vector_modulus(b)));
    return angle;
}

// 根据几何定理求夹角
double compute_steo_vector_with_inclination(Position a,Position b){
    SteoCord Sa = to_steocord(a);
    SteoCord Sb = to_steocord(b);
    SteoCord vec_1 = compute_steo_delta_vector(Sa,Sb);
    Position c = {a.latitude,b.longitude};
    SteoCord Sc = to_steocord(c);
    SteoCord vec_2 = compute_steo_delta_vector(Sa,Sc);
    double angle = compute_steo_vectors_angle(vec_1,vec_2);
    return angle;
}

// 确定最终的角度，以水平方向东北方向为正方向
double compute_steo_real_angle(Position a,Position b){
    if(b.longitude > a.longitude){
        double angle = compute_steo_vector_with_inclination(a,b);
        if(b.latitude>a.latitude)return angle;
        else return -angle;
    }
    else{
        double angle = compute_steo_vector_with_inclination(b,a);
        if(a.latitude>b.latitude)return angle;
        else return -angle;
    }
    return 0.0;
}

// 计算投影的长度
Distance compute_project_distance(Position a,Position b){
    Distance dist;
    double angle = compute_steo_real_angle(a,b);
    double delt_angle = angle - INCLINE_ANGLE;
    dist.true_dist = compute_curve_distance(a,b);
    dist.along_orbit = dist.true_dist*cos(to_radians(delt_angle));
    dist.vertical_orbit = dist.true_dist*sin(to_radians(delt_angle));
    return dist;
}

// 赤道平面的法向量为

// double compute_sphere_angle(struct Position x, struct Position y){
//     // 球面三角形的角度，待确定
//     double angle_betw = acos(sin(to_radians(x.latitude))*sin(to_radians(y.latitude))+
//                              cos(to_radians(x.latitude))*cos(to_radians(y.latitude))*cos(to_radians(y.longitude-x.longitude)));
//     return angle_betw;
// }

// struct Distance compute_all_distance(struct Position x, struct Position y){
//     struct Distance dist;
//     dist.true_dist = compute_haversine_distance(x,y);
//     double angle = compute_sphere_angle(x,y);
    
// }