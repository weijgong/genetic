/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-31 21:52:12
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2024-01-02 00:39:56
 * @FilePath: /root/genetic/genetic_v1/main.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 * 
 */
#include "coding.h"

using namespace std;

 
bool** gene_init_main_pop(){
    srand((unsigned)time(NULL));
    bool **pop = (bool **)malloc(MAX_TARGET_NUM * sizeof(bool *));
    int rand_num = 0;
    for(int i = 0;i < MAX_TARGET_NUM;i ++){
        rand_num=rand()%MainCoding_CASE;
        // cout<<rand_num<<endl;
        pop[i] = encode_main_code(rand_num);
        // for(int j = 0;j < MAIN_BIN_BITS;j ++){
        //     cout<<pop[i][j];
        // }
        // cout<<endl;
    }
    return pop;
}

bool* encode_main_code(int dec){
    int tmp = dec;
    bool* bin_arr = (bool*)malloc(sizeof(bool)*MAIN_BIN_BITS);
    for(int i = 0;i < MAIN_BIN_BITS;i ++){
        bin_arr[MAIN_BIN_BITS-1-i] = tmp%2;
        tmp/=2;
    }
    return bin_arr;

    // for(int i = 0;i < MAIN_BIN_BITS;i ++){
    //     cout<<bin_arr[i];
    // }
    // cout<<endl;

}

int decode_main_code(bool* bin_arr){
    int dec = 0;
    int mul = 1;

    for(int i = 0;i < MAIN_BIN_BITS;i ++){
        dec+=mul*bin_arr[MAIN_BIN_BITS-1-i];
        mul*=2;
    }
    // cout<<dec<<endl;
    return dec;
}

// 可以处理任意长度的十进制
void encode_bin_from_dec(int dec){
    int bin_bits = 0;
    int tmp = dec;
    while(tmp!=0){
        tmp/=2;
        bin_bits+=1;
    }
    cout<<bin_bits<<endl;
    bool* bin_arr = (bool*)malloc(sizeof(bool)*bin_bits);
    for(int i = 0;i < bin_bits;i ++){
        bin_arr[i] = 0;
    }
    tmp = dec;
    for(int i = 0;i < bin_bits;i ++){
        bin_arr[bin_bits-1-i] = tmp%2;
        tmp/=2;
    }
    for(int i = 0;i < bin_bits;i ++){
        cout<<bin_arr[i];
    }
    cout<<endl;
}

// int main(){
//     // bool* a = encode_main_code(3);
//     // int b = decode_main_code(a);
//     // cout<<b<<endl;
//     gene_init_main_pop();
// }