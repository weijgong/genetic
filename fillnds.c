/* Nond-domination based selection routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"

/* Routine to perform non-dominated sorting */
void fill_nondominated_sort (NSGA2Type *nsga2Params, population *mixed_pop, population *new_pop)
{
    int flag;
    int i, j;
    int end;
    int front_size;
    int archieve_size;
    int rank=1;
    list *pool;
    list *elite;
    list *temp1, *temp2;
    pool = (list *)malloc(sizeof(list));
    elite = (list *)malloc(sizeof(list));
    front_size = 0;
    archieve_size=0;
    pool->index = -1;
    pool->parent = NULL;
    pool->child = NULL;
    elite->index = -1;
    elite->parent = NULL;
    elite->child = NULL;
    temp1 = pool;
    // 将所有的个体索引存入pool中
    for (i=0; i<2*nsga2Params->popsize; i++)
    {
        insert (temp1,i);
        temp1 = temp1->child;
    }
    i=0;
    do
    {
        // ========================== 遍历pool中所有的个体 ==========================
        temp1 = pool->child;
        insert (elite, temp1->index);
        front_size = 1;
        temp2 = elite->child;
        temp1 = del (temp1);
        temp1 = temp1->child;
        do
        {   
            // 遍历其中的前沿集F_1
            temp2 = elite->child;
            if (temp1==NULL)break;            
            do
            {
                // ========================== 检测索引p与其余索引q之间的支配关系 ==========================
                end = 0;
                flag = check_dominance (nsga2Params, &(mixed_pop->ind[temp1->index]), &(mixed_pop->ind[temp2->index]));
                if (flag == 1)
                {
                    // p支配q,说明q不属于支配前沿,放入S_p中,检测下一个支配前沿
                    insert (pool, temp2->index);
                    temp2 = del (temp2);
                    front_size--;
                    temp2 = temp2->child;
                }
                if (flag == 0)
                {
                    // p与q互相不为支配解,检测下一个支配前沿
                    temp2 = temp2->child;
                }
                if (flag == -1)
                {
                    // p被q支配,找到了p的支配解,跳出循环
                    end = 1;
                }
            }
            while (end!=1 && temp2!=NULL);
            if (flag == 0 || flag == 1)
            {
                // 如果p支配所有的q或者成不支配的关系,将解p的索引放入精英集F_1中
                insert (elite, temp1->index);
                front_size++;
                temp1 = del (temp1);
            }
            temp1 = temp1->child;
        }
        while (temp1 != NULL);
        temp2 = elite->child;
        j=i;
        if ( (archieve_size+front_size) <= nsga2Params->popsize)
        {
            do
            {
                // ========================== 设定elite集中的各个元素rank=0 ==========================
                copy_ind (nsga2Params, &mixed_pop->ind[temp2->index], &new_pop->ind[i]);
                new_pop->ind[i].rank = rank;
                archieve_size+=1;
                temp2 = temp2->child;
                i+=1;
            }
            while (temp2 != NULL);
            assign_crowding_distance_indices (nsga2Params, new_pop, j, i-1);
            rank+=1;
        }
        else
        {
            crowding_fill (nsga2Params, mixed_pop, new_pop, i, front_size, elite);
            archieve_size = nsga2Params->popsize;
            for (j=i; j<nsga2Params->popsize; j++)
            {
                new_pop->ind[j].rank = rank;
            }
        }
        temp2 = elite->child;
        do
        {
            temp2 = del (temp2);
            temp2 = temp2->child;
        }
        while (elite->child !=NULL);
    }
    while (archieve_size < nsga2Params->popsize);
    // ========================== 释放两个list的空间 ==========================
    while (pool!=NULL)
    {
        temp1 = pool;
        pool = pool->child;
        free (temp1);
    }
    while (elite!=NULL)
    {
        temp1 = elite;
        elite = elite->child;
        free (temp1);
    }
    return;
}

/* Routine to fill a population with individuals in the decreasing order of crowding distance */
void crowding_fill (NSGA2Type *nsga2Params, population *mixed_pop, population *new_pop, int count, int front_size, list *elite)
{
    int *dist;
    list *temp;
    int i, j;
    assign_crowding_distance_list (nsga2Params,  mixed_pop, elite->child, front_size);
    dist = (int *)malloc(front_size*sizeof(int));
    temp = elite->child;
    for (j=0; j<front_size; j++)
    {
        dist[j] = temp->index;
        temp = temp->child;
    }
    quicksort_dist (mixed_pop, dist, front_size);
    for (i=count, j=front_size-1; i<nsga2Params->popsize; i++, j--)
    {
        copy_ind(nsga2Params,  &mixed_pop->ind[dist[j]], &new_pop->ind[i]);
    }
    free (dist);
    return;
}
