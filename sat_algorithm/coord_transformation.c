/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-04 10:27:34
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-04 10:36:54
 * @FilePath: /gongweijing/nsga2/sat_algorithm/coord_transformation.c
 * @Description: 经纬度转换为直角坐标系
 * 
 * Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
 */
#include "sat_algorithm.h"

SteoCord to_steocord(struct Position x){
    SteoCord cord;
    cord.x = EARTH_RADIUS_METER * cos(to_radians(x.latitude)) * cos(to_radians(x.longitude));
    cord.y = EARTH_RADIUS_METER * cos(to_radians(x.latitude)) * sin(to_radians(x.longitude));
    cord.z = EARTH_RADIUS_METER * sin(to_radians(x.latitude));
    return cord;
}