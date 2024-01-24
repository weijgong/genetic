/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-02 01:33:15
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2024-01-23 01:23:24
 * @FilePath: /gongweijing/genetic/sat_algorithm/read_access.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */

#include "sat_algorithm.h"
using namespace std;


string getExecuatePath(){
    char buff[FILENAME_MAX];
    getcwd( buff, FILENAME_MAX );
    string current_working_dir(buff);
    return current_working_dir;
}

string getExecuateParentPath(){
    string current_working_dir = getExecuatePath();
    size_t lastSlashPos = current_working_dir.find_last_of('/');

    if (lastSlashPos != std::string::npos) {
        // Extract the substring up to the last slash position
        std::string result = current_working_dir.substr(0, lastSlashPos+1);
        // std::cout << "Processed path: " << result << std::endl;
        return result;
    } else {
        // Handle the case when no slash is found
        std::cout << "Invalid path" << std::endl;
        return "";
    }
}

void parse_csv_line(char *line, struct AccessRecord *record) {
    char *token;
    token = strtok(line, ",");
    record->target_no = atoi(token);

    token = strtok(NULL, ",");
    strncpy(record->start_time, token, MAX_LINE_LENGTH);

    token = strtok(NULL, ",");
    strncpy(record->stop_time, token, MAX_LINE_LENGTH);
    
    token = strtok(NULL,",");
    record->duration = atoi(token);
}

void init_time_windows(){
        // 文件指针
    FILE *file;

    string root_path = getExecuateParentPath();
    // string data_path = "/sat_data/Satellite-Satellite1-Sensor-Sensor1-To-Target-Target1_Access.csv";
    // string data_path = "/sat_data/2024-1-22-Satellite-Satellite1-Sensor-Sensor1-To-Target-Target1_Access.csv";
    string data_path = "/sat_data/10-target-time-windows.csv";
    string filepath = root_path.append(data_path);
    // 文件路径
    const char *file_path = filepath.c_str();

    // 打开文件以供读取
    file = fopen(file_path, "r");

    // 检查文件是否成功打开
    if (file == NULL) {
        fprintf(stderr, "无法打开文件：%s\n", file_path);
        exit(1);
    }

    char line[MAX_LINE_LENGTH];

    int n = 0;
    // 读取CSV文件
    while (fgets(line, sizeof(line), file) != NULL) {
        struct AccessRecord record;
        parse_csv_line(line, &record);
        // printf("%d %s %s\n",record.target_no,record.start_time,record.stop_time);
        tw_list[n].target_no = record.target_no;
        tw_list[n].start_time=utc_to_tai(record.start_time);
        tw_list[n].stop_time=utc_to_tai(record.stop_time);
        tw_list[n].durations=tw_list[n].stop_time-tw_list[n].start_time;
        // printf("%ld %ld\n",utc_to_tai(record.start_time),utc_to_tai(record.stop_time));
        n++;
    }
    // 关闭文件
    fclose(file);
}

void extend_time_windows(){
    for(int i = 0;i<MAX_TARGET_NUM;i++){
        time_t cent_time = (tw_list[i].start_time+tw_list[i].stop_time)/2;
        tw_list[i].start_time = cent_time - tw_list[i].durations*5;
        tw_list[i].stop_time = cent_time + tw_list[i].durations*5;
        tw_list[i].durations *=2*5;
    }
}

// int main() {
//     init_time_windows();
//     return 0;  // 返回零表示程序执行成功
// }
