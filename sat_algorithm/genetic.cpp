#include "genetic.h"
using namespace std;

int GeneticAlgorithm::GetPopulationSize(){
    return this->populationSize;
}

Individual GeneticAlgorithm::reGenerateOfferingInfo(Individual ind,Sense_mode* SenseModeArray){
    for(int i = 0;i < MAX_TARGET_NUM;i ++){
        ind.genes[i].set_mode(SenseModeArray);
        ind.genes[i].exec_central_window();
    }
    ind = this->RepairUnfeasibleSolution(ind);
    ind.fitness = this->calculateFitness(ind);
    return ind;
}

GeneticAlgorithm::GeneticAlgorithm(int popSize, double mutationRate, double crossoverRate)
        : populationSize(popSize), mutationRate(mutationRate), crossoverRate(crossoverRate) {
        }


Individual GeneticAlgorithm::RepairUnfeasibleSolution(Individual ec){
    // b genetic.cpp:24
    vector<double> start_time_copy(MAX_TARGET_NUM);
    vector<int>    target_index(MAX_TARGET_NUM);
    for(int i=0;i<MAX_TARGET_NUM;i++){
        start_time_copy[i] = ec.genes[i].real_start_time;        
    }

    vector<pair<double, int>> vp;
    for (int i = 0; i < MAX_TARGET_NUM; ++i) {
        vp.push_back(make_pair(start_time_copy[i], i));
    }
    sort(vp.begin(), vp.end());
    for (int i = 0; i < vp.size(); i++) {
        target_index[i]=vp[i].second;
        // cout<<target_index[i]<<endl;
    }

    double occupied_timestamp = ec.genes[target_index[0]].real_finish_time;
    for(int i = 1;i < MAX_TARGET_NUM;i++){
        if(ec.genes[target_index[i]].pop_main_decode!=0 || ec.genes[target_index[i]].pop_main_decode!=3){
            if(occupied_timestamp > ec.genes[target_index[i]].real_start_time){
                ec.genes[target_index[i]].mode = -1;
                ec.genes[target_index[i]].pop_main_decode = 0;
                // 这个仅仅是对个体里面的一个基因进行了重新的位置编码
                ec.genes[target_index[i]].exec_central_window();
            }
            else{
                if(occupied_timestamp < ec.genes[target_index[i]].real_finish_time){
                    occupied_timestamp = ec.genes[target_index[i]].real_finish_time;
                }
            }
        }
    }

    vector<vector<double>>task_seting; 
    vector<vector<double>>mode3_task_seting; 
    double start = 0;
    int pos = 0;
    int free_tw_cnt=0;
    vector<vector<double>> interval_free;
    vector<double> vtw_s;
    vector<double> vtw_e;
    vector<double>tmps;
    for(int i=0;i<MAX_TARGET_NUM;i++){
        tmps.push_back(ec.genes[i].window.start_time);
        vtw_s.push_back(ec.genes[i].window.start_time);
        tmps.push_back(ec.genes[i].window.stop_time);
        vtw_e.push_back(ec.genes[i].window.stop_time);
        tmps.push_back(ec.genes[i].real_start_time);
        tmps.push_back(ec.genes[i].real_finish_time);
        tmps.push_back((double)ec.genes[i].pop_main_decode);
        tmps.push_back((double)i);
        task_seting.push_back(tmps);
        vector<double>().swap(tmps);
    }

    sort_ow_in_reg_Arrive(task_seting,task_seting.size());
    for(int i=0;i<MAX_TARGET_NUM;i++){
        if(task_seting[i][4]==0 || task_seting[i][4]==3)continue;
        if(task_seting[i][2]-start>=0.5){
            tmps.push_back(start);
            tmps.push_back(task_seting[i][2]);
            interval_free.push_back(tmps);
            vector<double>().swap(tmps);
            free_tw_cnt+=1;
        }
        start=task_seting[i][3];
    }
    tmps.push_back(start);
//  b genetic.cpp:94
    auto it_s = min_element(vtw_s.begin(), vtw_s.end());
    auto it_e = max_element(vtw_e.begin(), vtw_e.end());
    tmps.push_back((*it_e)-(*it_s));
    interval_free.push_back(tmps);
    vector<double>().swap(tmps);

    sort_in_reg_Arrive(task_seting,task_seting.size());
    vector<vector<double>> interval_free_resch = interval_free;
    vector<int>     unsch_task_index;                    // 存储未调度任务的下标
    for(int i=0;i<MAX_TARGET_NUM;i++){
        if(task_seting[i][4]==3)unsch_task_index.push_back(i);
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
            ec.genes[(int)task_seting[unsch_task_index[i]][5]].real_start_time  = task_seting[unsch_task_index[i]][2];
            ec.genes[(int)task_seting[unsch_task_index[i]][5]].real_finish_time = task_seting[unsch_task_index[i]][3];

            // task_seting[unsch_task_index[i]][4]=3;
            interval_free_resch[best_windows_index][0]+=0.5;
            if(interval_free_resch[best_windows_index][1]-interval_free_resch[best_windows_index][0]<=0.5){
                interval_free_resch.erase(interval_free_resch.begin()+best_windows_index);    
            }        
            unsch_task_index.erase(unsch_task_index.begin()+i);
        }
        else{
            ec.genes[(int)task_seting[unsch_task_index[i]][5]].pop_main_decode = 0;
            ec.genes[(int)task_seting[unsch_task_index[i]][5]].mode = -1;
        }
    }

    return ec;

    /* version 1:
    vector<double> start_time_copy(MAX_TARGET_NUM);
    vector<int>    target_index(MAX_TARGET_NUM);
    for(int i=0;i<MAX_TARGET_NUM;i++){
        start_time_copy[i] = ec.genes[i].real_start_time;        
    }

    vector<pair<double, int>> vp;
    for (int i = 0; i < MAX_TARGET_NUM; ++i) {
        vp.push_back(make_pair(start_time_copy[i], i));
    }
    sort(vp.begin(), vp.end());
    for (int i = 0; i < vp.size(); i++) {
        target_index[i]=vp[i].second;
        // cout<<target_index[i]<<endl;
    }

    double occupied_timestamp = ec.genes[target_index[0]].real_finish_time;
    for(int i = 1;i < MAX_TARGET_NUM;i++){
        if(ec.genes[target_index[i]].pop_main_decode!=0){
            if(occupied_timestamp > ec.genes[target_index[i]].real_start_time){
                ec.genes[target_index[i]].mode = -1;
                ec.genes[target_index[i]].pop_main_decode = 0;
                // 这个仅仅是对个体里面的一个基因进行了重新的位置编码
                ec.genes[target_index[i]].exec_central_window();
            }
            else{
                if(occupied_timestamp < ec.genes[target_index[i]].real_finish_time){
                    occupied_timestamp = ec.genes[target_index[i]].real_finish_time;
                }
            }
        }
    }
    return ec;
    */

}

// 获取随机基因
vector<EvaluationCode> GeneticAlgorithm::getRandomGene(Sense_mode *SenseModeArray) {
    vector<EvaluationCode>ec_simple(MAX_TARGET_NUM);
    for(int i = 0;i < MAX_TARGET_NUM; i++){
        // 种群初始化
        ec_simple[i].init_individual(SenseModeArray,i);
    }
    return ec_simple;
}

// 初始化种群
vector<Individual> GeneticAlgorithm::initializePopulation(Sense_mode *SenseModeArray) {
    vector<Individual> population;
    for (int i = 0; i < populationSize; ++i) {
        Individual ec_ind;
        ec_ind.genes = getRandomGene(SenseModeArray);
        ec_ind = this->RepairUnfeasibleSolution(ec_ind);
        ec_ind.fitness = this->calculateFitness(ec_ind);
        population.push_back(ec_ind);
    }
    return population;
}

void GeneticAlgorithm::nout_individual(Individual ec){
    cout<<"St\tFt\tfitness\tencode"<<endl;
    for(int i = 0;i < MAX_TARGET_NUM;i ++){
        printf("%.2f\t%.2f\t%.2f\t%d\n",
        ec.genes[i].real_start_time,ec.genes[i].real_finish_time,ec.fitness,ec.genes[i].pop_main_decode);
    }
}

void GeneticAlgorithm::abort_population(vector<Individual>Population){
    for(int i = 0;i < Population.size();i ++){
        cout<<"########### 第"<<i<<"个个体 ###########"<<endl;
        this->nout_individual(Population[i]);
    }
}

// 计算个体适应度
double GeneticAlgorithm::calculateFitness(Individual ec) {
    double fitness = 0;
    int observed_target_num=0;
    for(int i = 0; i < MAX_TARGET_NUM; i ++){
        if(ec.genes[i].pop_main_decode==0){            
        }
        else if(ec.genes[i].pop_main_decode==1){
            fitness+=10;
            observed_target_num+=1;
        }
        else if(ec.genes[i].pop_main_decode==2){
            fitness+=5;
            observed_target_num+=1;
        }
        else if(ec.genes[i].pop_main_decode==3){
            fitness+=1;
            observed_target_num+=1;
        }
    }
    fitness += observed_target_num;
    return fitness;
}

void GeneticAlgorithm::assginFitness(vector<Individual>& Population){
    for(int i = 0;i < Population.size();i ++){
        Population[i].fitness = this->calculateFitness(Population[i]);
    }
}

// 从pool中找出最好的一个个体
Individual GeneticAlgorithm::Tournament(vector<Individual> TournamentCandidate){
    Individual best = TournamentCandidate[0];
    for(int i = 1; i < TournamentCandidate.size();i ++){
        if(best.fitness < TournamentCandidate[i].fitness){
            best = TournamentCandidate[i];
        }
    }
    return best;
}

// 随机取出不重复的pool_number个个体
vector<Individual> GeneticAlgorithm::SampleTournament(vector<Individual> population,int pool_number){
    vector<int> _index(this->populationSize,0);
    vector<Individual> TnmSet;
    for(int i = 0;i < this->populationSize;i ++){
        _index[i] = i;
    }
    random_device rd;
    mt19937 g(rd());
    shuffle(_index.begin(),_index.end(),g);
    // cout<<"seq\treal_index"<<endl;
    for(int i = 0;i < pool_number;i ++){
        // printf("%d\t%d\n",i,_index[i]);
        TnmSet.push_back(population[_index[i]]);
    }    
    return TnmSet;
}

// 锦标赛算法进行选择
vector<Individual> GeneticAlgorithm::TournamentSelection(vector<Individual> population,
                                                         int select_number,
                                                         int tournament_pool_size
                                                         ){
    if(tournament_pool_size > this->populationSize){
        printf("Error: 锦标赛选择尺寸过大,自动修正为种群大小");
        tournament_pool_size = this->populationSize;
    }                                                            
    vector<Individual>select_pop;
    vector<Individual>tournament_pool_pop;
    for(int i = 0;i < select_number;i ++){
        tournament_pool_pop = this->SampleTournament(population,tournament_pool_size);
        select_pop.push_back(this->Tournament(tournament_pool_pop));
    }
    return select_pop;
}

// 交叉
vector<Individual> GeneticAlgorithm::crossover(vector<Individual> parent) {
    double randmask = (float) rand() / RAND_MAX;
    vector<int> mask;
    vector<Individual>child;
    child = parent;

    for(int i = 0;i < MAX_TARGET_NUM;i ++){
        randmask = (float) rand() / RAND_MAX;
        int randi;
        if(randmask>0.5)randi=1;
        else randi = 0;
        mask.push_back(randi);
    }

    for(int i = 0;i < MAX_TARGET_NUM;i ++){
        if(mask[i]==1){
            child[0].genes[i] = parent[1].genes[i];
            child[1].genes[i] = parent[0].genes[i];
        }
    }

    for(int i=0;i<child.size();i++){
        child[i].fitness = this->calculateFitness(child[i]);
    }

    return child;

}

// 变异
Individual GeneticAlgorithm::mutate(Individual ind) {
    for (int i = 0;i < MAX_TARGET_NUM;i ++) {
        if ( (float) rand() / RAND_MAX < this->mutationRate) {
            // 整数编码中的任意一个都算是符合要求的
            ind.genes[i].pop_main_decode = rand()%MainCoding_CASE;
        }
    }
    return ind;
}
