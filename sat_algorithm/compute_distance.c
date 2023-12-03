/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-03 21:16:52
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-03 21:41:30
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