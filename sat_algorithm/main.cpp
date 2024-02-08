/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-02 01:33:21
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2024-01-23 14:26:11
 * @FilePath: /gongweijing/genetic/sat_algorithm/main.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
// I.E.
#include "sat_algorithm.h"
#include "coding.h"
#include "genetic.h"

using namespace std;
using namespace std::chrono;
void nout_individual_gene(Individual ind);


struct TimeWindow tw_list [MAX_TARGET_NUM];
struct Position   pos_list[MAX_TARGET_NUM];
struct SteoCord   cor_list[MAX_TARGET_NUM];



bool compare_by_idx(vector<double> a,vector<double> b,int idx){
    if(a[idx]>=b[idx]){
        return true;
    }
    else return false;
}

void sort_in_reg_Arrive(vector<vector<double>> &arr,int size_){
    quick_sort_by_arrive(arr,0,size_-1);
}

void sort_ow_in_reg_Arrive(vector<vector<double>> &arr,int size_){
    quick_sort_by_ow_arrive(arr,0,size_-1);
}

void sort_target_no(vector<vector<double>> &arr,int size_){
    quick_sort_by_target_no(arr,0,size_-1);
}

void quick_sort_by_target_no(vector<vector<double>> &arr, int start, int end) {
    if (start >= end)
        return;
    vector<double> mid = arr[end];
    int left = start, right = end - 1;
    while (left < right) {
        while (!compare_by_idx(arr[left],mid,5) && left < right)
            left++;
        while (compare_by_idx(arr[right],mid,5) && left < right)
            right--;
        swap(arr[left], arr[right]);
    }
    if (compare_by_idx(arr[left],arr[end],5))
        swap(arr[left], arr[end]);
    else
        left++;
    quick_sort_by_target_no(arr, start, left - 1);
    quick_sort_by_target_no(arr, left + 1, end);
}

void quick_sort_by_ow_arrive(vector<vector<double>> &arr, int start, int end) {
    if (start >= end)
        return;
    vector<double> mid = arr[end];
    int left = start, right = end - 1;
    while (left < right) {
        while (!compare_by_idx(arr[left],mid,2) && left < right)
            left++;
        while (compare_by_idx(arr[right],mid,2) && left < right)
            right--;
        swap(arr[left], arr[right]);
    }
    if (compare_by_idx(arr[left],arr[end],2))
        swap(arr[left], arr[end]);
    else
        left++;
    quick_sort_by_ow_arrive(arr, start, left - 1);
    quick_sort_by_ow_arrive(arr, left + 1, end);
}

void quick_sort_by_arrive(vector<vector<double>> &arr, int start, int end) {
    if (start >= end)
        return;
    vector<double> mid = arr[end];
    int left = start, right = end - 1;
    while (left < right) {
        while (!compare_by_idx(arr[left],mid,0) && left < right)
            left++;
        while (compare_by_idx(arr[right],mid,0) && left < right)
            right--;
        swap(arr[left], arr[right]);
    }
    if (compare_by_idx(arr[left],arr[end],0))
        swap(arr[left], arr[end]);
    else
        left++;
    quick_sort_by_arrive(arr, start, left - 1);
    quick_sort_by_arrive(arr, left + 1, end);
}

void nout_individual_gene(Individual ind){
    cout<<"ind code:\t";
    for(int i = 0; i < MAX_TARGET_NUM;i ++){
        cout<<ind.genes[i].pop_main_decode;
    }cout<<endl;
}

int main(){
    
    auto start_time = high_resolution_clock::now();
    
    srand((unsigned)time(NULL));
    Sense_mode *SenseModeArray = NULL;

    vector<int>           _mode(MAX_TARGET_NUM);
    vector<vector<double>> vtws(MAX_TARGET_NUM,vector<double>(2));
    vector<vector<double>> ows (MAX_TARGET_NUM,vector<double>(2));

    int numMode = 0;
    read_sense_mode(&SenseModeArray,&numMode);
    init_time_windows();

    extend_time_windows();
    time_t earliest_time_start = (time_t)INFINITY;
    time_t slowes_time_stop = 0;
    for(int i = 0;i<MAX_TARGET_NUM;i++){
        if(earliest_time_start>tw_list[i].start_time){
            earliest_time_start = tw_list[i].start_time;            
        }
        if(slowes_time_stop<tw_list[i].stop_time){
            slowes_time_stop = tw_list[i].stop_time;
        }
    }

    // cout<<"最早的时间窗口开始时间为："<<earliest_time_start<<endl;
    // cout<<"最晚的时间窗口结束时间为："<<slowes_time_stop<<endl;

    // 主要作用是进行约束条件和编码后生成开始窗口
    for(int i = 0;i<MAX_TARGET_NUM;i++){
        tw_list[i].start_time -= earliest_time_start;
        tw_list[i].stop_time  -= earliest_time_start;

        // printf("Start:%ld,end:%ld,duration:%ld\n",tw_list[i].start_time,tw_list[i].stop_time,tw_list[i].durations);
    }


    // *********************************** Evolve Process ***********************************
    int populationSize = 100;
    double mutation_ratio = 0.1;
    double crossover_ratio = 0.9;
    int generation_number = 200;

    GeneticAlgorithm algo(populationSize,mutation_ratio,crossover_ratio);
    vector<Individual> population = algo.initializePopulation(SenseModeArray);
    Individual best_ind = algo.Tournament(population);
    algo.nout_individual(best_ind);

    int poolNumber = 2;
    
    for(int epoch = 0;epoch < generation_number; epoch ++){
        vector<Individual> newPopulation(populationSize);
        Individual cur_best_ind;
        Individual best_child;
        for(int i = 0;i < algo.GetPopulationSize();i ++){
            vector<Individual> child(poolNumber);
            vector<Individual> parent = algo.TournamentSelection(population,2,6);

            child = algo.crossover(parent);

            child[0] = algo.mutate(child[0]);
            child[1] = algo.mutate(child[1]);

            child[0] = algo.reGenerateOfferingInfo(child[0],SenseModeArray);
            child[1] = algo.reGenerateOfferingInfo(child[1],SenseModeArray);
            
            best_child=algo.Tournament(child);
            newPopulation[i] = best_child;
        }
      
        cur_best_ind = algo.Tournament(newPopulation);

        if(best_ind.fitness < cur_best_ind.fitness){
            // 每当更新最优的个体的时候就abort一次
            printf("第%d次",epoch);
            cout<<"出现更优的个体"<<&cur_best_ind<<",更新\n";
            best_ind = cur_best_ind;
            algo.nout_individual(best_ind);
        }
        
        population = newPopulation;
        newPopulation.clear();

    }
    
    cout<<"最好的个体"<<&best_ind<<"："<<endl;
    algo.nout_individual(best_ind);

    // *********************************** Evolve Process End ***********************************
    FILE *vtw_ow_fileP;
    vtw_ow_fileP = fopen("/home/gwj/genetic/sat_algorithm/time_windows_data.dat", "w");
    if (vtw_ow_fileP == nullptr) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    for(int i=0;i<MAX_TARGET_NUM;i++){
        vtws[i][0] = best_ind.genes[i].window.start_time;
        vtws[i][1] = best_ind.genes[i].window.stop_time;
        ows [i][0] = best_ind.genes[i].real_start_time;
        ows [i][1] = best_ind.genes[i].real_finish_time;
        _mode  [i] = best_ind.genes[i].pop_main_decode;
        fprintf(vtw_ow_fileP, "%f,%f,%f,%f,%d\n",vtws[i][0],vtws[i][1],ows [i][0],ows [i][1],_mode  [i]);
    }
    
    fclose(vtw_ow_fileP);

    
    FILE *free_twP;
    free_twP = fopen("/home/gwj/genetic/sat_algorithm/free_time_windows_data.dat", "w");
    
    if (free_twP == nullptr) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // plot_target_windows("all_access_timewindows.png",earliest_time_start,slowes_time_stop);
    // plot_individual_with_name("best_pop_in_simple_genetic.png",best_ind,earliest_time_start,slowes_time_stop);
    // plot_individual_on_oneline_with_name("best_pop_in_simple_genetic_oneline.png",best_ind,earliest_time_start,slowes_time_stop);


    // *********************************** Re Schedule Process Start ***********************************
    vector<vector<double>>task_seting;                                   // (vtw_s,vtw_e,ow_s,ow_e,mode,target_no)
    
    double start = 0;
    int pos = 0;
    int free_tw_cnt=0;
    vector<vector<double>> interval_free;
    vector<double>tmps;
    for(int i=0;i<MAX_TARGET_NUM;i++){
        tmps.push_back(vtws[i][0]);
        tmps.push_back(vtws[i][1]);
        tmps.push_back(ows[i][0]);
        tmps.push_back(ows[i][1]);
        tmps.push_back((double)_mode[i]);
        tmps.push_back((double)i);
        task_seting.push_back(tmps);
        vector<double>().swap(tmps);
    }

    // vector<vector<double>>task_seting_by_no = task_seting;   

    sort_ow_in_reg_Arrive(task_seting,task_seting.size());
    for(int i=0;i<MAX_TARGET_NUM;i++){
        if(task_seting[i][4]==0)continue;
        if(task_seting[i][2]-start>=0.5){
            tmps.push_back(start);
            tmps.push_back(task_seting[i][2]);
            interval_free.push_back(tmps);
            vector<double>().swap(tmps);
            fprintf(free_twP, "%f,%f\n",interval_free[free_tw_cnt][0],interval_free[free_tw_cnt][1]);
            free_tw_cnt+=1;
        }
        start=task_seting[i][3];
    }
    tmps.push_back(start);
    tmps.push_back(slowes_time_stop-earliest_time_start);
    interval_free.push_back(tmps);
    vector<double>().swap(tmps);
    fprintf(free_twP, "%f,%f\n",interval_free[free_tw_cnt][0],interval_free[free_tw_cnt][1]);


    // 三种排序方法：到达时间优先、截止时间优先、等待时间优先
    // 按照vtws开始时间进行排序
    // Regulation AI:
    sort_in_reg_Arrive(task_seting,task_seting.size());
    vector<vector<double>> interval_free_resch = interval_free;
    vector<int>     unsch_task_index;                    // 存储未调度任务的下标
    vector<int>       sch_task_index;                    // 存储调度任务的下标
    for(int i=0;i<MAX_TARGET_NUM;i++){
        if(task_seting[i][4]==0)unsch_task_index.push_back(i);
        else sch_task_index.push_back(i);
    }

    for(int i=unsch_task_index.size()-1;i>=0;i--){
        double cur_s,cur_e;
        double best_wci=0,cur_wci=0;
        double best_windows_index = -1;
        cur_s = task_seting[unsch_task_index[i]][0];
        cur_e = task_seting[unsch_task_index[i]][1];
        for(int j=0;j<interval_free_resch.size();j++){
            cur_wci = min(cur_e,interval_free_resch[j][1]) - max(cur_s,interval_free_resch[j][0]);
            if(cur_wci >= 0.5){
                // 满足条件的时间窗口属于正常的
                // incident_windows[i].push_back(j);
                if(cur_wci>=best_wci){
                    best_wci = cur_wci;
                    best_windows_index = j;
                }
            }
        }
        if(best_windows_index!=-1){
            task_seting[unsch_task_index[i]][2]=max(interval_free_resch[best_windows_index][0],
                                                    task_seting[unsch_task_index[i]][2]);
            task_seting[unsch_task_index[i]][3]=task_seting[unsch_task_index[i]][2]+0.5;
            task_seting[unsch_task_index[i]][4]=3;
            interval_free_resch[best_windows_index][0]+=0.5;
            if(interval_free_resch[best_windows_index][1]-interval_free_resch[best_windows_index][0]<=0.5){
                interval_free_resch.erase(interval_free_resch.begin()+best_windows_index);    
            }        
            unsch_task_index.erase(unsch_task_index.begin()+i);
        }
    }
    
    FILE *rescheduled_vtw_ow_fileP;
    rescheduled_vtw_ow_fileP = fopen("/home/gwj/genetic/sat_algorithm/rescheduled_time_windows_data.dat", "w");

    if (rescheduled_vtw_ow_fileP == nullptr) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }
    sort_target_no(task_seting,task_seting.size());
    for(int i=0;i<MAX_TARGET_NUM;i++){
        for(int j=0;j<6;j++){
            printf("%f\t",task_seting[i][j]);
        }
        fprintf(rescheduled_vtw_ow_fileP, "%f,%f,%f,%f,%f\n",
        task_seting[i][0],
        task_seting[i][1],
        task_seting[i][2],
        task_seting[i][3],
        task_seting[i][4]);
        printf("\n");
    }

    // *********************************** Re Schedule Process End ***********************************
    auto end_time = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end_time - start_time);
 
    cout << "Time taken by function: "
         << duration.count() << " microseconds" << endl;
    fclose(vtw_ow_fileP);
    free(SenseModeArray);
    return 0;
// */
}

/*
    vector<double> wci_task_i_windows_j;
    for(int i=0;i<incident_windows.size();i++){
        double cur_s,cur_e;
        cur_s = task_seting[unsch_task_index[i]][0];
        cur_e = task_seting[unsch_task_index[i]][1];
        // cout<<"incident_windows[i] size is "<<incident_windows[i].size()<<"\t";
        // cout<<"target no is "<<task_seting[unsch_task_index[i]][5]<<":\t";
        for(int j=0;j<incident_windows[i].size();j++){
            cout<<incident_windows[i][j]<<"\t";
            cout<<"wci->"<<min(cur_e,interval_free[incident_windows[i][j]][1]) - max(cur_s,interval_free[incident_windows[i][j]][0])<<"\t";
            wci_task_i_windows_j.push_back(min(cur_e,interval_free[j][1]) - max(cur_s,interval_free[j][0]));
        }
        cout<<"\n";
    }
*/

/*
    for(int i=0;i<interval_free.size();i++){
        printf("%f,%f\n",interval_free[i][0],interval_free[i][1]);
    }

    vector<vector<vector<double>>> incident_windows(unsch_task_index.size(),
    vector<vector<double>>(interval_free.size(),vector<double>(2)));
    printf("incident windows size=%d,unsch_task_index windows size=%d\n",
        incident_windows.size(),unsch_task_index.size());
    for(int i=0;i<incident_windows.size();i++)
        for(int j=0;j<interval_free.size();j++)
            for(int k=0;k<2;k++)
                incident_windows[i][j][k] = interval_free[j][k];
    for(int i=0;i<incident_windows.size();i++){
        for(int j=0;j<incident_windows[0].size();j++){
            printf("%f~%f\t",incident_windows[i][j][0],incident_windows[i][j][1]);
        }
        printf("\n");
    }


    for(int i=0;i<incident_windows.size();i++){
        sort(deleted_windows[i].begin(), deleted_windows[i].end(), greater<size_t>());
        for (size_t index : deleted_windows[i]) {
            if (index < incident_windows[i].size()) {
                incident_windows[i].erase(incident_windows[i].begin() + index);
            }
        }
    }

    for(int i=0;i<incident_windows.size();i++){
        printf("%f:\t", task_seting[unsch_task_index[i]][5]);
        for(int j=0;j<incident_windows[0].size();j++){
            printf("%f~%f\t",incident_windows[i][j][0],incident_windows[i][j][1]);
        }
        printf("\n");
    }
*/ 

/*
    plot_target_windows("all_access_timewindows.png",earliest_time_start,slowes_time_stop);
    plot_individual_with_name("best_pop_in_simple_genetic.png",best_ind,earliest_time_start,slowes_time_stop);
    plot_individual_on_oneline_with_name("best_pop_in_simple_genetic_oneline.png",best_ind,earliest_time_start,slowes_time_stop);
*/

/*
    for(int i=0;i<MAX_TARGET_NUM;i++){
        for(int j=0;j<6;j++){
            printf("%f\t",task_seting[i][j]);
        }
        printf("\n");
    }
*/
    
/*
    for (int i = 0; i < numMode; i++) {
        printf("Mode %d:\n", i + 1);
        printf("Imaging Duration: %.2lf seconds\n", SenseModeArray[i].imagingDuration);
        printf("Resolution: %.2lf meters\n", SenseModeArray[i].resolution);
        printf("Swath Width: %.2lf kilometers\n", SenseModeArray[i].swathWidth);
        printf("Is Slant View: %d\n", SenseModeArray[i].isSlantView);
        printf("Slant Angle: %.2lf degrees\n", SenseModeArray[i].slantAngle);
        printf("\n");
    }  

    for(int i = 0;i<MAX_TARGET_NUM;i++){
        printf("Target no %d:\n", tw_list[i].target_no);
        printf("Start time windows: %ld \n", tw_list[i].start_time);
        printf("End time windows: %ld \n", tw_list[i].stop_time);
        printf("Time windows duration: %ld \n", tw_list[i].durations);
        printf("\n");
    }

    printf("Target no %d:\n", tw_list[i].target_no);
    printf("Start time windows: %ld \n", tw_list[i].start_time);
    printf("End time windows: %ld \n", tw_list[i].stop_time);
    printf("Time windows duration: %ld \n", tw_list[i].durations);
    printf("\n");

    init_position();
    for(int i = 0;i<MAX_TARGET_NUM;i++){
        printf("Target Latitude: %f \n", pos_list[i].latitude);
        printf("Target Longitude: %f \n", pos_list[i].longitude);
        printf("\n");
    }
*/

/*
    cout<<"初始最优个体"<<&best_ind<<"\n";
    algo.nout_individual(best_ind);
    
    algo.nout_individual(population[0]);
    plot_individual_with_name("test_plot.png",population[0],earliest_time_start,slowes_time_stop);
    plot_individual(population[0],earliest_time_start,slowes_time_stop);
    algo.abort_population(population);
*/
    
/*
    cout<<"pop_num\tfitness\tencoding\n";
    for(int i = 0;i < algo.GetPopulationSize();i ++){
        printf("%d\t%f\t",i,population[i].fitness);
        for(int mg = 0;mg < MAX_TARGET_NUM;mg ++){
            printf("%d ",population[i].genes[mg].pop_main_decode);
        }cout<<endl;
    }
*/

/*
    cout<<"最好的个体："<<endl;
    algo.nout_individual(best_ind);
    plot_individual_with_name("best_pop_in_simple_genetic.png",best_ind,earliest_time_start,slowes_time_stop);
    
    Individual binds = algo.Tournament(population);
    algo.nout_individual(binds);

    // 求解最好的个体就对所有的进行锦标赛，找出其中最好的个体
    Individual best_ind = algo.Tournament(population);
    cout<<"最好的个体："<<endl;
    algo.nout_individual(best_ind);
    plot_individual_with_name("test_plot.png",best_ind,earliest_time_start,slowes_time_stop);
*/

/*
    // ************** 测试代码👇 **************
    vector<Individual> newPopulation;
    vector<Individual> parent = algo.TournamentSelection(population,2,6);
    // 输出生成的两个父代,便于进行交叉重组
    algo.abort_population(parent);

    vector<Individual> child = algo.crossover(parent);
    // cout<<"子代个体交叉后:"<<endl;
    // algo.abort_population(child);

    algo.mutate(child[0]);
    algo.mutate(child[1]);

    // cout<<"子代个体变异后:"<<endl;
    // algo.abort_population(child);

    algo.reGenerateOfferingInfo(child[0]);
    algo.reGenerateOfferingInfo(child[1]);

    // cout<<"子代个体修复后:"<<endl;
    algo.abort_population(child);

    algo.assginFitness(child);

    // cout<<"适应度函数重新计算后:"<<endl;
    // algo.abort_population(child);

    Individual best_child = algo.Tournament(child);
    // algo.abort_population(child);
    algo.nout_individual(best_child);
    plot_individual_with_name("test_plot.png",best_child,earliest_time_start,slowes_time_stop);
    // ************** 测试代码👆 **************
*/
    


        /*
            // 问题大概率在下面这部分
            gdb main  b main.cpp:185 r 
            display parent
            display child
            parent被浅拷贝进了child,由于parent本身设定的是一个不断被删除生成的数据,所以随着地址的内容不断变化,进而影响了剩余全部的数据;
            在crossover进行深拷贝之后，截断了child跟parent指向同一个EvaluationCode，后续也就不会出现覆盖写的问题
            p algo.abort_population (parent)
            display nout_individual_gene(child[0])
        */
    /*
        cout<<"循环外\n";
        
        printf("%d\t%f\t",epoch,cur_best_ind.fitness);
        for(int mg = 0;mg < MAX_TARGET_NUM;mg ++){
            printf("%d ",cur_best_ind.genes[mg].pop_main_decode);
        }
        cout<<endl;
    */
        // cout<<"当前种群最好的个体:"<<&cur_best_ind<<endl;
        // algo.nout_individual(cur_best_ind);
/*
    cout<<"Epoch\tBestFit\t\tEncoding\n";
    
        cout<<"pop_num\tfitness\t\tencoding\n";
        
            // b main.cpp:191 display best_child  display newPopulation display child
            // display nout_individual_gene(*best_child) 
            // 问题查出：上一次放入population的变量被child变量覆盖了
            // 解决方案，直接一次性生成n个个体，然后每次赋值

            // 确定是否为repair的bug
            algo.abort_population(child);

            cout<<"第"<<i<<"轮的子代对个体地址为："<<&(child[0])<<" "<<&(child[1])<<endl;
            cout<<"第"<<i<<"轮的子代对最优个体地址为："<<&best_child<<endl;
            cout<<"存储到Pupulation的地址为"<<&(newPopulation[i])<<endl;

            cout<<"第"<<i<<"轮的子代对\n";
            algo.nout_individual(*best_child);

            cout<<"第"<<i<<"轮的子代对\n";
            algo.abort_population(*child);

            Individual* best_child = new Individual;
            *best_child = algo.Tournament(*child);
            newPopulation.push_back(*best_child);
            algo.nout_individual(*best_child);

            if(newPopulation.size()!=0){
                cout<<"第i轮循环变化为："<<endl;
                algo.nout_individual(newPopulation[i-1]);
            }

            cout<<"before cross\n";
            algo.abort_population(parent);

            cout<<"before repair\n";
            algo.abort_population(child);
            cout<<"第"<<i<<"轮的子代对\n";
            algo.nout_individual(child[0]);
            for(int mg = 0;mg < MAX_TARGET_NUM; mg ++){
                child[0].genes[mg].exec_central_window();
            }

            cout<<"after repair\n";
            algo.abort_population(child);

            cout<<"第i-1轮送入数据为："<<endl;
            algo.nout_individual(newPopulation[i]);


            algo.nout_individual(best_child);
            printf("%d\t%f\t",i,best_child.fitness);
            for(int mg = 0;mg < MAX_TARGET_NUM;mg ++){
                printf("%d ",best_child.genes[mg].pop_main_decode);
            }cout<<endl;


            if(i == populationSize-1){
                cout<<"循环内:\n";
                cur_best_ind = algo.Tournament(newPopulation);
                algo.nout_individual(cur_best_ind);
            }
*/
          