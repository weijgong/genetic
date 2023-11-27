/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-11-25 09:39:53
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-11-26 09:45:43
 * @FilePath: /gongweijing/nsga2/auxiliary.c
 * @Description: 
 * 
 * Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
 */
/* Some utility functions (not part of the algorithm) */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"

/* Function to return the maximum of two variables */
double maximum (double a, double b)
{
    if (a>b)
    {
        return(a);
    }
    return (b);
}

/* Function to return the minimum of two variables */
double minimum (double a, double b)
{
    if (a<b)
    {
        return (a);
    }
    return (b);
}
