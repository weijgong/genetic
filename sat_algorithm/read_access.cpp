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

void nout_struct_tm(struct tm tm_){
    printf("Year: %d\n", tm_.tm_year + 1900); // 年份
    printf("Month: %d\n", tm_.tm_mon + 1);     // 月份
    printf("Day: %d\n", tm_.tm_mday);          // 日
    printf("Hour: %d\n", tm_.tm_hour);         // 小时
    printf("Minute: %d\n", tm_.tm_min);        // 分钟
    printf("Second: %d\n", tm_.tm_sec);        // 秒
}

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
    string root_path = getExecuateParentPath();
    // string data_path = "/sat_data/Satellite-Satellite1-Sensor-Sensor1-To-Target-Target1_Access.csv";
    // string data_path = "/sat_data/2024-1-22-Satellite-Satellite1-Sensor-Sensor1-To-Target-Target1_Access.csv";
    string data_path = "sat_data/10-target-time-windows.csv";
    string filepath = root_path.append(data_path);

    vector<string> start_s;
    vector<string> stop_s;
    vector<string> duration_s;
    vector<string> target_no_s;

    ifstream inputFile(filepath); // 打开文件
    if (inputFile.is_open()) { // 检查文件是否成功打开
        string line;

        int n = 0;
        int pos = 0;
        while (std::getline(inputFile, line)) {
            pos = line.find(',');
            target_no_s.push_back(line.substr(0,pos));
            line = line.substr(pos+1,line.size());
            pos = line.find(',');
            start_s.push_back(line.substr(0,pos));
            line = line.substr(pos+1,line.size());
            pos = line.find(',');
            stop_s.push_back(line.substr(0,pos));        
            line = line.substr(pos+1,line.size());
            pos = line.find(',');
            duration_s.push_back(line.substr(0,pos));
            
            n+=1;
        }
        inputFile.close();
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
    for(int i=0;i<MAX_TARGET_NUM;i++){
        tw_list[i].target_no=stoi(target_no_s[i]);
        tw_list[i].start_time=tm_to_seconds(start_s[i]);
        tw_list[i].stop_time =tm_to_seconds(stop_s[i]);
        tw_list[i].durations=stoi(duration_s[i]);

        // printf("info:%ld,%ld,%ld,%ld\n",
        // tw_list[i].target_no,
        // tw_list[i].start_time,
        // tw_list[i].stop_time,
        // tw_list[i].durations);

        // printf("%s,%s,%s,%s\n",target_no_s[i].data(),start_s[i].data(),stop_s[i].data(),duration_s[i].data());

        // cout<<tm_to_seconds(start_s[i])<<"\n";
        // cout<<tm_to_seconds(stop_s[i])<<"\n";
    }
}

void extend_time_windows(){
    for(int i = 0;i<MAX_TARGET_NUM;i++){
        time_t cent_time = (tw_list[i].start_time+tw_list[i].stop_time)/2;
        tw_list[i].start_time = cent_time - tw_list[i].durations*5;
        tw_list[i].stop_time = cent_time + tw_list[i].durations*5;
        tw_list[i].durations *=2*5;
        // printf("%d,%d,%d\n",tw_list[i].start_time,tw_list[i].stop_time,tw_list[i].durations);
    }
}

// int main() {
//     init_time_windows();
//     return 0;  // 返回零表示程序执行成功
// }
