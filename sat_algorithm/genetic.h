/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2024-01-02 00:21:19
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2024-01-02 00:25:02
 * @FilePath: /root/genetic/sat_algorithm/genetic.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%A
 */
#ifndef GENETIC_H
#define GENETIC_H
#include "coding.h"

struct Individual {
    EvaluationCode* ec_simple;
    double fitness;
};

#endif