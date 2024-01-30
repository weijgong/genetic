/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-11-25 09:39:53
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-11-29 09:42:54
 * @FilePath: /gongweijing/nsga2/report.c
 * @Description:   
 * 
 * Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
 */
/* Routines for storing population data into files */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"

/**
 * @description: 依次将目标函数数值、约束条件数值、实数变量值、二进制变量值、约束冲突、解的等级、拥挤距离存入.out文件
 * @param {NSGA2Type} *nsga2Params
 * @param {population} *pop
 * @param {FILE} *fpt
 * @return {*}
 */
void report_pop (NSGA2Type *nsga2Params,  population *pop, FILE *fpt)
{
    int i, j, k;
    for (i=0; i<nsga2Params->popsize; i++)
    {
        for (j=0; j<nsga2Params->nobj; j++)
        {
            fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
        }
        if (nsga2Params->ncon!=0)
        {
            for (j=0; j<nsga2Params->ncon; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
            }
        }
        if (nsga2Params->nreal!=0)
        {
            for (j=0; j<nsga2Params->nreal; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].xreal[j]);
            }
        }
        if (nsga2Params->nbin!=0)
        {
            for (j=0; j<nsga2Params->nbin; j++)
            {
                for (k=0; k<nsga2Params->nbits[j]; k++)
                {
                    fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                }
            }
        }
        fprintf(fpt,"%e\t",pop->ind[i].constr_violation);
        fprintf(fpt,"%d\t",pop->ind[i].rank);
        fprintf(fpt,"%e\n",pop->ind[i].crowd_dist);
    }
    return;
}

/**
 * @description: 输出所有的可行解,可行解的定义为：pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1
 * @param {NSGA2Type} *nsga2Params
 * @param {population} *pop
 * @param {FILE} *fpt
 * @return {*}
 */
void report_feasible (NSGA2Type *nsga2Params,  population *pop, FILE *fpt)
{
    int i, j, k;
    for (i=0; i<nsga2Params->popsize; i++)
    {
        if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1)
        {
            for (j=0; j<nsga2Params->nobj; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
            }
            if (nsga2Params->ncon!=0)
            {
                for (j=0; j<nsga2Params->ncon; j++)
                {
                    fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
                }
            }
            if (nsga2Params->nreal!=0)
            {
                for (j=0; j<nsga2Params->nreal; j++)
                {
                    fprintf(fpt,"%e\t",pop->ind[i].xreal[j]);
                }
            }
            if (nsga2Params->nbin!=0)
            {
                for (j=0; j<nsga2Params->nbin; j++)
                {
                    for (k=0; k<nsga2Params->nbits[j]; k++)
                    {
                        fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                    }
                }
            }
            fprintf(fpt,"%e\t",pop->ind[i].constr_violation);
            fprintf(fpt,"%d\t",pop->ind[i].rank);
            fprintf(fpt,"%e\n",pop->ind[i].crowd_dist);
        }
    }
    return;
}
