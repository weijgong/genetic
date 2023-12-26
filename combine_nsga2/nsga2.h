/*
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-20 19:32:37
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-20 19:32:59
 * @FilePath: /root/genetic/combine_nsga2/nsga2.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <unistd.h>

# define INF 1.0e14
# define EPS 1.0e-14
# define eul  2.71828182845905
# define pi 3.14159265358979
# define GNUPLOT_COMMAND "gnuplot -persist"


typedef struct
{
    int rank;
    double constr_violation;
    double *xreal;
    int **gene;
    double *xbin;
    double *obj;
    double *constr;
    double crowd_dist;
}
individual;

typedef struct
{
    individual *ind;
}
population;

typedef struct lists
{
    int index;
    struct lists *parent;
    struct lists *child;
} list;


typedef struct NSGA2Type
{
    double seed;
    int nreal;
    int nbin;
    int nobj;
    int ncon;
    int popsize;
    double pcross_real;
    double pcross_bin;
    double pmut_real;
    double pmut_bin;
    double eta_c;
    double eta_m;
    int ngen;
    int nbinmut;
    int nrealmut;
    int nbincross;
    int nrealcross;
    int *nbits;
    double *min_realvar;
    double *max_realvar;
    double *min_binvar;
    double *max_binvar;
    int bitlength;
    int choice;
    int obj1;
    int obj2;
    int obj3;
    int angle1;
    int angle2;
} NSGA2Type;

// Global 
extern FILE *fpt1;
extern FILE *fpt2;
extern FILE *fpt3;
extern FILE *fpt4;
extern FILE *fpt5;
extern FILE *gp;
extern population *parent_pop;
extern population *child_pop;
extern population *mixed_pop;
