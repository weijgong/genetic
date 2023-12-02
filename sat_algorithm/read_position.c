/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-03 00:17:15
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-03 00:51:42
 * @FilePath: /root/genetic/sat_algorithm/read_position.c
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "sat_algorithm.h"

void init_position(){
    FILE *file;
    // 文件路径
    const char *file_path = "/root/genetic/sat_data/Target1_LLA_Position.csv";
    // 打开文件以供读取
    file = fopen(file_path, "r");
    // 检查文件是否成功打开
    if (file == NULL) {
        fprintf(stderr, "无法打开文件：%s\n", file_path);
        return 1;  // 返回非零值表示程序出错
    }
    char line[MAX_LINE_LENGTH];
    int n = 0;
    // 读取CSV文件
    while (fgets(line, sizeof(line), file) != NULL) {
        struct AccessRecord record;
        char *token;
        token = strtok(line, ",");
        pos_list[n].latitude = atoi(token);

        token = strtok(line, ",");
        pos_list[n++].longitude = atoi(token);       
    }
}
int main() {
    init_position();
    return 0;  // 返回零表示程序执行成功
}