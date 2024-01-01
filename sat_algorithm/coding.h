/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2024-01-01 15:23:38
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2024-01-01 17:34:33
 * @FilePath: /root/genetic/genetic_v1/coding.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef CODING_H
#define CODING_H

#include "sat_algorithm.h"

#include <iostream>
#include <ctime>
#include <cstdlib>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>


#define MAIN_BIN_BITS 2
#define MAX_TARGET_NUM 5

// 三种模式加上一种什么模式都不选择的
#define MainCoding_CASE 4


class EvaluationCode{
    public:
        // parameter
        int target_no;
        double real_start_time;
        double real_finish_time;
        int exec_satellite_index;
        struct TimeWindow window;
        bool* pop_main;
        int * pop_sub;        
        int fitness;
        int pop_main_decode;

        int mode;
        double observe_time;

        // function
        void exec_central_window();
        void nout_time_proc();
        void set_mode(Sense_mode *SenseModeArray);
};

// 生成初始种群
bool** gene_init_main_pop();
// 十进制转二进制
bool* encode_main_code(int dec);
// 二进制转十进制
int decode_main_code(bool* bin_arr);
// 可以处理任意长度的十进制
void encode_bin_from_dec(int dec);


#endif