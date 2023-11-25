/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-11-25 09:39:53
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-11-25 15:32:29
 * @FilePath: /gongweijing/nsga2/list.c
 * @Description: list insert/Delete operate.
 * 
 * Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
 */
/* A custom doubly linked list implemenation */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"

/**
 * @description: Insert the list node with index x after the node node
 * @param {list} *node
 * @param {int} x
 * @return {*}
 */
void insert (list *node, int x)
{
    list *temp;
    if (node==NULL)
    {
        printf("\n Error!! asked to enter after a NULL pointer, hence exiting \n");
        exit(1);
    }
    temp = (list *)malloc(sizeof(list));
    temp->index = x;
    temp->child = node->child;
    temp->parent = node;
    if (node->child != NULL)
    {
        node->child->parent = temp;
    }
    node->child = temp;
    return;
}

/* Delete the node NODE from the list */
list* del (list *node)
{
    list *temp;
    if (node==NULL)
    {
        printf("\n Error!! asked to delete a NULL pointer, hence exiting \n");
        exit(1);
    }
    temp = node->parent;
    temp->child = node->child;
    if (temp->child!=NULL)
    {
        temp->child->parent = temp;
    }
    free (node);
    return (temp);
}
