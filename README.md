This is the Readme file for NSGA-II code.


About this version
==================
Basically this version is a *Refactored* version of the original code in order to make the code structure more portable. Several changes have been made:

1. The **NSGA2Type** type has been defined to carry all the necessary parameters of the NSGA2 algorithm. This will reduce the number of _external_ variables and reduce the possible conflict when it uses within another routine.
2. Three function has been extracted from the main() function, previously in **nsga2r.c**. 
	1. ReadParamters: Read the parameters through input file or command line
	2. InitNSGA2: Initialize output files, allocate memory, and perform the first generation
	3. NSGA2: Run the generation, save the results and free the allocated memory.
	4. PrintNSGA2Parameters: Print the input parameters.
3. Two extra void pointers have been provided and passed around to the function that has to call `evaluate` function. These are helpful when you want to integrate the algorithm to your code (e.g. simulator). Using `void *inp` and `void *out` pointer you could pass your simulator parameters around without changing the structure of NSGA2, as long as your objective function knows how to handle the `inp` and `out` to produce the objective values.

About the Algorithm
===================
NSGA-II: Non-dominated Sorting Genetic Algorithm - II

Please refer to the following paper for details about the algorithm:

- Authors: Dr. Kalyanmoy Deb, Sameer Agrawal, Amrit Pratap, T Meyarivan
- Paper Title: A Fast and Elitist multi-objective Genetic Algorithm: NSGA-II
- Journal: IEEE Transactions on Evolutionary Computation (IEEE-TEC)
- Year: 2002
- Volume: 6
- Number: 2
- Pages: 182-197


NOTE
====

This archive contains routines for plotting the objective data real-time using gnuplot. The code has been written for posix compliant operating systems and uses sta
ndard piping method provided by GNU C library. The routines should work on any unix and unix like OS having gnuplot installed and which are posix compliant.



How to compile and run the program
==================================
Makefile has been provided for compiling the program on linux (and unix-like) systems. Edit the Makefile to suit your need. By default, provided Makefile attempts to compile and link all the existing source files into one single executable.

Name of the executable produced is: nsga2r

To run the program type: `./nsga2r random_seed` or `./nsga2r random_seed < input_data/inp_file.in`

Here random_seed is a real number in (0,1) which is used as a seed for random number generator, and "inp_file.in" is the file that stores all the input parameters


About the output files
======================

| File | Description |
|:----------|:------------|
| initial_pop.out | This file contains all the information about initial population. |
| final_pop.out | This file contains the data of final population. |
| all_pop.out | This file contains the data of populations at all generations. |
| best_pop.out | This file contains the best solutions obtained at the end of simulation run. |
| params.out | This file contains the information about input parameters as read by the program. |



About the input parameters
==========================

| Parameter | Description |
|:----------|:------------|
| popsize | This variable stores the population size (a multiple of 4) |
| ngen | Number of generations |
| nobj | Number of objectives |
| ncon | Number of constraints |
| nreal | Number of real variables |
| min_realvar[i] | minimum value of i^{th} real variable |
| max_realvar[i] | maximum value of i^{th} real variable |
| pcross_real | probability of crossover of real variable |
| pmut_real | probability of mutation of real variable |
| eta_c | distribution index for real variable SBX crossover |
| eta_m | distribution index for real variable polynomial mutation |
| nbin | number of binary variables |
| nbits[i] | number of bits for i^{th} binary variable |
| min_binvar[i] | minimum value of i^{th} binary variable |
| max_binvar[i] | maximum value of i^{th} binary variable |
| pcross_bin | probability of crossover for binary variable |
| pmut_bin | probability of mutation for binary variable |
| choice | option to display the data realtime using gnuplot |
| obj1, obj2, obj3 | index of objectives to be shown on x, y and z axes respectively |
| angle1, angle2 | polar and azimuthal angle required for location of eye |



定义测试问题
=========================
编辑源文件 problemdef.c 以定义测试问题。我们提供了一些示例问题（Deb 博士的书 - 使用进化算法的多目标优化中的 24 个测试问题）作为示例，以指导您如何定义自己的目标和约束函数。您还可以根据需要将其他源文件与代码链接。在实现目标函数时，请考虑以下几点。

1. 代码是为最小化目标而编写的（最小f_i）。如果要最大化函数，可以使用函数值的负数作为目标值。
2. 如果一个解决方案不违反任何约束，则该解决方案被认为是可行的。如果解必须可行，约束函数的计算结果应大于或等于零 （g_j >= 0）。约束的负值意味着它被违反。
3. 如果存在多个约束条件，建议（尽管不是强制性的）通过重新制定约束值或将其除以正非零常数来规范约束值。

About the files
===============

| Files | Descriptions |
|:-----------|:------------|
| nsga2.h | 包含全局变量和函数声明的头文件 |
| rand.h | 包含随机数生成器的变量和函数声明的头文件 |
| allocate.c | 内存分配和释放例程 |
| auxiliary.c | 辅助例程（不属于算法的一部分） |
| crossover.c | 实数和二进制交叉的例程 |
| crowddist.c | 拥挤距离分配例程 |
| decode.c | 解码二进制变量的例程 |
| display.c | 使用 gnuplot 实时显示数据的例程 |
| dominance.c | 执行非支配检查的例程 |
| eval.c | 用于评估约束冲突的例程 |
| fillnds.c | 基于非支配排序的选择 |
| initialize.c | 对填充成员执行随机初始化的例程 |
| list.c | 自定义双向链表实现 |
| merge.c | 将两个种群合并为一个较大种群的例程 |
| mutation.c | 实数突变和二元突变的例程 |
| nsga2r.c | 主要功能的实现和NSGA-II框架 |
| problemdef.c | 测试问题定义 |
| rand.c | 随机数生成器相关例程 |
| rank.c | 排名分配例程 |
| report.c | 在文件中写入人口信息的例程 |
| sort.c | 随机快速排序实现 |
| tourselect.c | 锦标赛选择程序 |


Please feel free to send questions/comments/doubts/suggestions/bugs
etc. to deb@iitk.ac.in

Dr. Kalyanmoy Deb
14th June 2005
http://www.iitk.ac.in/kangal/

