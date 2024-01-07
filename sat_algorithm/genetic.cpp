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
}

// 获取随机基因
EvaluationCode* GeneticAlgorithm::getRandomGene(Sense_mode *SenseModeArray) {
    EvaluationCode* ec_simple = (EvaluationCode*)malloc(sizeof(EvaluationCode)*MAX_TARGET_NUM);
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
    for(int i = 0; i < MAX_TARGET_NUM; i ++){
        if(ec.genes[i].pop_main_decode==0){            
        }
        else if(ec.genes[i].pop_main_decode==1){
            fitness+=10;
        }
        else if(ec.genes[i].pop_main_decode==2){
            fitness+=5;
        }
        else if(ec.genes[i].pop_main_decode==3){
            fitness+=1;
        }
    }
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
    vector<Individual>child(2);

    // cout<<"mask:\t";
    for(int i = 0;i < MAX_TARGET_NUM;i ++){
        randmask = (float) rand() / RAND_MAX;
        int randi;
        if(randmask>0.5)randi=1;
        else randi = 0;
        mask.push_back(randi);
        // cout<<randi<<"\t";
    }
    // cout<<endl;

    Individual parent1 = parent[0];
    Individual parent2 = parent[1];
    // 赋值之前要先初始化空间出来，当然也能直接采用push_back加进去两个父代便于交叉
    child[0] = parent1;
    child[1] = parent2;
    
    for(int i = 0;i < MAX_TARGET_NUM;i ++){
        if(mask[i]==1){
            // 注意进行交叉的时候只对于每个编码进行交叉，其余的基础信息不需要进行交叉处理
            child[0].genes[i].pop_main_decode = parent2.genes[i].pop_main_decode;
            child[1].genes[i].pop_main_decode = parent1.genes[i].pop_main_decode;
        }
        else{
            // 好像初始化的时候对应的父代基因自动复制进去了，不需要进行赋值
        }
    }

    return child;

// #### 不确定是不是这样写，可能不是，先暂定不是
    // if(randmask<this->crossoverRate){
    //     vector<int> mask;
    //     vector<Individual>child(2);

    //     for(int i = 0;i < MAX_TARGET_NUM;i ++){
    //         randmask = (float) rand() / RAND_MAX;
    //         mask.push_back(randmask);
    //     }
    //     Individual parent1 = parent[0];
    //     Individual parent2 = parent[1];
    //     // 赋值之前要先初始化空间出来，当然也能直接采用push_back加进去两个父代便于交叉
    //     child[0] = parent1;
    //     child[1] = parent2;
        
    //     for(int i = 0;i < MAX_TARGET_NUM;i ++){
    //         if(mask[i]==1){
    //             // 注意进行交叉的时候只对于每个编码进行交叉，其余的基础信息不需要进行交叉处理
    //             child[0].genes[i].pop_main_decode = parent2.genes[i].pop_main_decode;
    //             child[1].genes[i].pop_main_decode = parent1.genes[i].pop_main_decode;
    //         }
    //         else{
    //             // 好像初始化的时候对应的父代基因自动复制进去了，不需要进行赋值
    //         }
    //     }

    //     return child;
    // }
    // else{
    //     return parent;
    // }
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

// 遗传算法主循环
// Individual evolve(int generations) {
//     std::vector<Individual> population = initializePopulation();

//     for (int generation = 0; generation < generations; ++generation) {
//         for (Individual& individual : population) {
//             calculateFitness(individual);
//         }

//         std::vector<Individual> newPopulation;

//         for (int i = 0; i < populationSize; ++i) {
//             Individual parent1 = selectParent(population);
//             Individual parent2 = selectParent(population);

//             Individual child = crossover(parent1, parent2);
//             mutate(child);

//             newPopulation.push_back(child);
//         }

//         population = newPopulation;
//     }

//     // 返回最终种群中适应度最高的个体
//     return *std::max_element(population.begin(), population.end(),
//                                 [](const Individual& ind1, const Individual& ind2) {
//                                     return ind1.fitness < ind2.fitness;
//                                 });
// }


// // 获取指定范围内的随机数
// double getRandomNumber(double min, double max) {
//     return min + static_cast<double>(rand()) / RAND_MAX * (max - min);
// }
