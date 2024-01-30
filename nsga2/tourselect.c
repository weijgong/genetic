/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-11-25 09:39:53
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-11-29 10:43:32
 * @FilePath: /gongweijing/nsga2/tourselect.c
 * @Description: 
 * 
 * Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
 */
/* Tournamenet Selections routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"

/**
 * @description: 用于锦标赛选择的例程，它通过执行锦标赛选择和交叉，从old_pop创建new_pop。
 * @param {NSGA2Type} *nsga2Params
 * @param {population} *old_pop
 * @param {population} *new_pop
 * @return {*}
 */
void selection (NSGA2Type *nsga2Params,  population *old_pop, population *new_pop)
{
    int *a1, *a2;
    int temp;
    int i;
    int rand;
    individual *parent1, *parent2;
    a1 = (int *)malloc(nsga2Params->popsize*sizeof(int));
    a2 = (int *)malloc(nsga2Params->popsize*sizeof(int));
    for (i=0; i<nsga2Params->popsize; i++)
    {
        a1[i] = a2[i] = i;
    }
    // 将数组索引进行随机打乱,便于后续进行随机tournament comparision.
    for (i=0; i<nsga2Params->popsize; i++)
    {
        rand = rnd (i, nsga2Params->popsize-1);
        temp = a1[rand];
        a1[rand] = a1[i];
        a1[i] = temp;
        rand = rnd (i, nsga2Params->popsize-1);
        temp = a2[rand];
        a2[rand] = a2[i];
        a2[i] = temp;
    }
    // 相邻两个个体比较，选出其中更优的一个进行交叉,每次生成4个新个体
    for (i=0; i<nsga2Params->popsize; i+=4)
    {
        parent1 = tournament (nsga2Params, &old_pop->ind[a1[i]], &old_pop->ind[a1[i+1]]);
        parent2 = tournament (nsga2Params, &old_pop->ind[a1[i+2]], &old_pop->ind[a1[i+3]]);
        crossover (nsga2Params, parent1, parent2, &new_pop->ind[i], &new_pop->ind[i+1]);
        parent1 = tournament (nsga2Params, &old_pop->ind[a2[i]], &old_pop->ind[a2[i+1]]);
        parent2 = tournament (nsga2Params, &old_pop->ind[a2[i+2]], &old_pop->ind[a2[i+3]]);
        crossover (nsga2Params, parent1, parent2, &new_pop->ind[i+2], &new_pop->ind[i+3]);
    }
    free (a1);
    free (a2);
    return;
}

/* Routine for binary tournament */
individual* tournament (NSGA2Type *nsga2Params,  individual *ind1, individual *ind2)
{
    int flag;
    flag = check_dominance (nsga2Params, ind1, ind2);
    if (flag==1)
    {
        return (ind1);
    }
    if (flag==-1)
    {
        return (ind2);
    }
    if (ind1->crowd_dist > ind2->crowd_dist)
    {
        return(ind1);
    }
    if (ind2->crowd_dist > ind1->crowd_dist)
    {
        return(ind2);
    }
    if ((randomperc()) <= 0.5)
    {
        return(ind1);
    }
    else
    {
        return(ind2);
    }
}
