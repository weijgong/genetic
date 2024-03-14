<!-- translated on 2023.11.27 -->
This is the Readme file for NSGA-II code.

pip3 install -i https://pypi.tuna.tsinghua.edu.cn/simple pandas
pip3 install -i https://pypi.tuna.tsinghua.edu.cn/simple scikit-learn
pip3 install -i https://pypi.tuna.tsinghua.edu.cn/simple scipy
About this version
==================
基本上，这个版本是原始代码的“重构”版本，以使代码结构更具可移植性。已经做了一些改变:
1. 已经定义了**NSGA2Type**类型来携带NSGA2算法所需的所有参数。这将减少_external_变量的数量，并减少在另一个例程中使用它时可能发生的冲突。
2. 在**nsga2r.c**中，从main()函数中提取了三个函数。
	1. readparameters:通过输入文件或命令行读取参数
	2. InitNSGA2:初始化输出文件，分配内存，执行第一代
	3. NSGA2:运行生成，保存结果并释放分配的内存。
	4. PrintNSGA2Parameters:打印输入参数。

3.提供了两个额外的void指针，并传递给必须调用`evaluate`函数的函数。当你想要将算法集成到你的代码(例如模拟器)时，这些是有用的。使用`void *inp`和`void *out`指针，你可以在不改变NSGA2结构的情况下传递模拟器参数，只要你的目标函数知道如何处理`inp` 和 `out`来产生目标值。

About the Algorithm
===================
NSGA-II:非支配排序遗传算法II (Non-dominated Sorting Genetic Algorithm - II)

有关算法的详细内容，请参考以下文件:
- 作者:Kalyanmoy Deb博士，Sameer Agrawal, Amrit Pratap, T Meyarivan
- 论文题目:一种快速的精英多目标遗传算法:NSGA-II
- 期刊:IEEE进化计算汇刊(IEEE- tec)
- 年份:2002年
- 容量:6
- 数量:2
- 页:182-197

NOTE
====
此存档包含使用gnuplot实时绘制客观数据的例程。该代码是为posix兼容的操作系统编写的，并使用GNU C库提供的标准管道方法。这些例程应该可以在任何安装了gnuplot并且兼容posix的unix和类似unix的操作系统上工作。


How to compile and run the program
==================================
Makefile提供了在linux(和类unix)系统上编译程序的功能。编辑Makefile以满足您的需要。默认情况下，提供Makefile尝试编译并将所有现有源文件链接到一个可执行文件中。

生成的可执行文件名称为:nsga2r

要运行程序，输入: `./nsga2r random_seed` or `./nsga2r random_seed < input_data/inp_file.in`

这里random_seed是(0,1)中的一个实数，它被用作随机数生成器的种子，而`inp_file.inp`是存储所有输入参数的文件


About the output files
======================
| File | Description |
|:----------|:------------|
| initial_pop.out | 这个文件包含所有关于初始种群的信息。|
| final_pop.out | 该文件包含最终种群的数据。|
| all_pop.out | 该文件包含所有世代种群的数据。|
| best_pop.out | 该文件包含模拟运行结束时获得的最佳解决方案。|
| params.out | 该文件包含程序读取的有关输入参数的信息。|


About the input parameters
==========================
| Parameter | Description |
|:----------|:------------|
| popsize |该变量存储总体大小(4的倍数)|
| n | 世代数 |
| nobj |目标数|
| ncon |约束数|
| nreal |实变量数|
| min_realvar[i] | $i^{th}$实变量的最小值|
| max_realvar[i] | $i^{th}$实变量的最大值|
| pcross_real | 实变量交叉的概率 |
| pmut_real | 实数变量的变异概率 |
| eta_c | 实变量SBX交叉的分布指标 |
| eta_m | 实数变量多项式突变的分布指标 |
| nbin |二进制变量数|
| nbits[i] | $i^{th}$二进制变量的位数|
| min_binvar[i] | $i^{th}$二进制变量的最小值|
| max_binvar[i] |第 $i^{th}$个二进制变量的最大值|
| pcross_bin |二元变量交叉概率
| pmut_bin |二进制变量的变异概率|
| choice | 选项显示数据实时使用gnuplot |
| obj1、obj2、obj3 |分别显示在x、y、z轴上的物镜索引|
| angle1, angle2 |眼睛定位所需的极角和方位角|

 


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

