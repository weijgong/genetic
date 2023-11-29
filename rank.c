/* Rank assignment routine */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"

/* Function to assign rank and crowding distance to a population of size pop_size*/
void assign_rank_and_crowding_distance (NSGA2Type *nsga2Params, population *new_pop)
{
    int flag;
    int i;
    int end;
    int front_size;
    int rank=1;
    list *orig;
    list *cur;
    list *temp1, *temp2;
    orig = (list *)malloc(sizeof(list));
    cur = (list *)malloc(sizeof(list));
    front_size = 0;
    orig->index = -1;
    orig->parent = NULL;
    orig->child = NULL;
    cur->index = -1;
    cur->parent = NULL;
    cur->child = NULL;
    temp1 = orig;
    // 初始化orig种群集的双向链表
    for (i=0; i<nsga2Params->popsize; i++)
    {
        insert (temp1,i);
        temp1 = temp1->child;
    }
    do
    {
        // 如果到达了orig最后一个元素,那么将对应元素的crowd_distance赋值为INF
        if (orig->child->child == NULL)
        {
            new_pop->ind[orig->child->index].rank = rank;
            new_pop->ind[orig->child->index].crowd_dist = INF;
            break;
        }
        // temp1为origin原本的第一个种群，将其插入cur中，并删除origin原本的第一个种群
        temp1 = orig->child;
        insert (cur, temp1->index);
        front_size = 1;
        temp2 = cur->child;
        temp1 = del (temp1);
        temp1 = temp1->child;
        // 第一层循环：origin不断遍历到child，直到到达空集结束；
        do
        {
            temp2 = cur->child;
            // 第二层循环：1)cur不断遍历到child，直到到达空集结束；2)origin对应节点被cur对应节点支配结束；
            do
            {
                end = 0;
                flag = check_dominance (nsga2Params, &(new_pop->ind[temp1->index]), &(new_pop->ind[temp2->index]));
                if (flag == 1)
                {
                    // 将cur内的节点放入origin，并从cur中将其删除，前沿数-1.
                    insert (orig, temp2->index);
                    temp2 = del (temp2);
                    front_size--;
                    temp2 = temp2->child;
                }
                if (flag == 0)
                {
                    temp2 = temp2->child;
                }
                if (flag == -1)
                {
                    end = 1;
                }
            }
            while (end!=1 && temp2!=NULL);
            if (flag == 0 || flag == 1)
            {
                // 将origin中不被cur上个结点支配的解放入cur中,作为候选前沿
                insert (cur, temp1->index);
                front_size++;
                temp1 = del (temp1);
            }
            temp1 = temp1->child;
        }
        while (temp1 != NULL);
        temp2 = cur->child;
        do
        {
            new_pop->ind[temp2->index].rank = rank;
            temp2 = temp2->child;
        }
        while (temp2 != NULL);
        assign_crowding_distance_list (nsga2Params, new_pop, cur->child, front_size);
        temp2 = cur->child;
        do
        {
            temp2 = del (temp2);
            temp2 = temp2->child;
        }
        while (cur->child !=NULL);
        rank+=1;
    }
    while (orig->child!=NULL);
    free (orig);
    free (cur);
    return;
}
