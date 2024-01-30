<!--
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-11-29 15:23:26
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-20 09:26:37
 * @FilePath: /gongweijing/nsga2/Project_Design.md
 * @Description: 记录项目主要框架
 * 
 * Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
-->
# 项目估计
感觉这个项目虽然功能实现的比较经典，但是由于接口有点复杂，每次配置的时候需要重新设定新的in文件+problem部分不能直接选择，
后续根据job shop schedule问题实现相关的项目找个遗传算法来进行实现，目前时间已经到了12.11，所以需要加快工作进度，防止
工期延误。

# 测试问题
## 多目标优化问题描述：生产计划问题

### 决策变量：
- $(x_1$): 生产产品1的数量（整数变量）
- $(x_2$): 生产产品2的数量（整数变量）

### 目标函数：
1. **最大化利润：**
   $[ Z_1 = - (20x_1 + 30x_2) ]$

2. **最小化产品数量：**
   $[ Z_2 = x_1 + x_2 ]$

### 约束条件：
$[ 2x_1 + 3x_2 \leq 60 ]$\
$[ 4x_1 + 2x_2 \leq 80 ]$\
$[ x_1 \leq 25 ]$\
$[ x_2 \leq 30 ]$\
$[ x_1, x_2 \geq 0 ] $（非负约束）


# Design 1
## object function:
1. 成像效果:
$$
f=\sum w_1(Solar\_angle - w_a) +  w_2(w_c- Coudy)
$$
其中，$Solar\_angle、Coudy $  分别表示成像任务所在位置的太阳高度角和云量等级;$w_a、w_c$分别表示太阳高度角下限和云量等级上限; $w_1,w_2$用来统一量纲,无实际意义。

2. 资源消耗:
$$
f= sum -\sum imageChange \times za/vl
$$
式中，imageChange为执行成像任务所需的侧视角度;za为侧视机动中单位时间内消耗的能量;vl为侧视机动的角速率。
为了保证多目标优化方向的一致性,资源消耗评价值采用了总能量与消耗能量之差。

3. 成像目标等级：
$$
value =\sum vertex\_value +dl\_qz;\\
vertex\_value=\begin{cases}
vertex\_qz,& d=0\\
\frac{vertex\_qz}{vertex\_qz^-}& d\ne 0
\end{cases}\\
vertex\_qz^-=\sum vertex\_qz^{\sim}\times vertex\_hc
$$
表达式中包括了成像任务的重要性权值和数传任务的重要性权值，$dl\_qz$为数传任务的重要性权值。

## constrain condition:
1. 成像任务时间窗口约束为$d(imagePosition_i,imagePosition_j)<6$。约束中$imagePosition_i$表示成像任务目标的地理位置,d表示距离函数。
    
    在经纬度标下计算两个任务目标之间的经纬度距离,当距离小于6时,由于卫星的飞行速度和成像开关机所需时间等原因,导致两个任务不能同时执行,只能任选其一。
    
2. 各地面站与卫星的通信仰角约束为$a_{ij}>5$。其中au 是指第i个卫星与第j个地面站的通信仰角。

    一般来说,当卫星与地面站之间的通信仰角小于5时,卫星与地面站就不能进行数据传输。

3. 通信星与成像星的通信距离约束为$d(Track\_S_{a_i},Track\_S_{b_i})<300km$。其中,$Track\_S_{a_i}$表示成像星在i时刻的经纬度坐标;$Track\_S_{b_i}$表示通信星在i时刻的经纬度坐标;
    
    当双星距离大于300km,认为不满足通信条件。

4. 可见光遥感器的光照条件约束为:$Solar\_angle_i<w_{a_i}$。约束公式中的$Solar\_angle_i$表示第i个成像任务目标在当前时刻的太阳高度角,$w_{a_i}$表示执行该任务需要满足的太阳高度角下限。

太阳高度角的计算公式如下:
$$
\delta= 23.45sin(\frac{2\pi d}{365})\\
Solar\_angle=sin\phi sin\delta+cos\phi cos\delta cos\omega
$$
d表示某年中某一天的日期序号;$\phi$表示当前成像任务目标的地理纬度;$\delta$表示太阳赤纬;$\omega$表示太阳时角。

5. 成像任务的云量等级约束为$Coudy_i < w_{c_i}$, 当不满足此约束条件,表明成像效果没有到达预期。$Coudy_i$表示成像任务目标所在位置上空的云量等级; $w_{c_i}$表示执行该成像任务所需的云量等级条件。

## math model
具有星间通信能力的成像卫星星座任务规划与传统多星任务规划有所不同,各任务之间的耦合性更强。不仅需要考虑成像任务之间的冲突,还需要考虑数传任务之间以及两种任务之间的冲突。成像星在其固定轨道上飞行,当其遥感器的覆盖范围出现了成像目标,就可能执行成像任务,当满足与地面站的通信条件,则可能执行数传任务。而当通信星的飞行轨道与成像星的轨道在极地区会合时,满足星问数传条件则可能执行星间数传任务。

设计星座中两颗星的轨道高度为600km,轨道类型为太阳同步轨道。通信星与成像星的星间数传距离约为200～300km,相应的星间相对相位约为1.7°~2.5°。

星间相对相位公式为：
$$\alpha =2arcsin(\frac{d}{h})$$
其中,$\alpha$为星间相对相位,d为星间相对距离(m),h为轨道高度(m)。

## STK simulate 1
### target(name longitude latitude)
1. 北京 39.9118 116.379
2. 海口 20.02   110.35
3. 喀什 39.48   76
4. 哈尔滨 45.76  126.62
### satellite(name 轨道高度(km) 轨道倾角(deg) 升交点(deg) 真近点角(deg))
1. 成像星 600 97.6402 0   350
2. 通信星 600 97.6402 100 1
### sensor and other informations
- 成像星携带的成像载荷为Siple Conic
- 坐标系统为J2000
- 仿真开始时间为:1th Jul 2007 12:00:00000,结束时间为:3th Jul 2007 12:00:00000,仿真步长为60sec.根据这 48hours 的仿真,可以得到每个采样时刻成像星和通信星在惯性坐标系下的位置信息和星下点轨迹信息。

## STK simulate 2
### satellite config
name: YAOGAN-35 04A
source: https://celestrak.org/satcat/search.php

### point target
```c
name,latitude,longitude
Point_Target1,6.443,110.537
Point_Target2,57.222,-20.094
Point_Target3,50.78,-5.025
Point_Target4,83.353,161.73
Point_Target5,-34.05,42.56
Point_Target6,79.04,176.275
Point_Target7,57.935,139.833
Point_Target8,-15.303,-79.731
Point_Target9,-28.67,128.083
Point_Target10,25.246,135.527
Point_Target11,-5.681,177.644
Point_Target12,-28.476,148.493
Point_Target13,-33.143,-143.135
Point_Target14,-5.18,-37.027
```

# param config
## 遗传算法
- 种群的初始数目可以设定为：200
- 种群的迭代次数可以设定为：300
- 交叉概率：0.9
- 变异概率：0.1
### 编码方式

### 约束配置

# Question:
- what is J2000?
- 


# 多目标
## 观测目标
- 尽可能完成多个目标的观测：
$$
f_1(x)=min\ -x_n \ 最大化观测数目\\
f_2(x)=min\ -\sum_i (h_i\times 10+m_i\times3+l_i\times 1) \ 观测收益\\
f_3(x)=min\ \sum (pic_1\times 10+pic_2 \times 3+pic_3 \times 1)\ 成像质量,对应三种成像模式；也跟成像角度有关，后续添加。
\\
$$

## 约束条件
### 时间窗口相关
$$
ta_i表示观测开始时间\\
dt_i表示成像时长,与成像模式相关\\
\\
\begin{cases}
T_{start}\le ta_i& 最早观测开始时间不早于可见窗口开始时刻T_{start}\\
ta_i\le T_{end}-dt_i& 最晚观测开始时间不晚于可见窗口结束时刻T_{end}-成像时长\\
ta_{i+1}>ta_i& 下一次成像时刻晚于当前成像时刻\\
\end{cases}
$$
### 成像区域相关约束
$$
pos_i,pos_j表示两观测目标的经纬度\\
dist表示计算两个目标的球面距离\\
amp表示成像幅宽\\
\begin{cases}
dist_{pos_i,pos_j}\le amp& \\
dist_{pos_i,pos_j}\le \frac{amp}{2}& 表明两个目标可以一次性被观测\\

\end{cases}
$$
## mode

## 输出
target(观测的目标)|mode(选择的模式)|start time(观测开始时间MKT)
|--|--|--|
|8|0|3441241|

## 仿真数据
- 时间窗口长度:40~60秒
- 窗口特点：1）相互间隔；2）相互重叠；
- 成像目标: 
    1. 个数：5个
    1. 重要性（0/1/2）; 
    2. 经度，纬度；小于一半幅宽，大于一半幅宽。
    3. 卫星轨道倾角35deg，计算投影的时候计算垂直轨道方向。

## SAR：
1. 入射角范围：12~52


## 思路设计
### 变量：
1. 是否观测：$x_i$， 表示第 $i$ 个目标是否被观测（取值为0或1）。
3. 观测目标可用时间窗口： $T_i$，表示的是第$i$个目标可以被观测到的时间窗口,具体有成像开始、结束时间、成像时间间隔。
2. 成像模式：$m_i$， 表示第 $i$ 个目标采用的观测模式，如果不进行观测默认赋值为-1。
5. 观测模式：$z_i$， 表示第 $i$ 个目标采用的观测模式（0表示正视，1表示侧视），如果不进行观测默认赋值为-1。
4. 成像幅宽：$amp_i$，表示第 $i$ 个目标进行观测的时候采用的具体成像幅宽，如果不进行观测默认赋值为-1。
6. 成像所需时间：$dt_i$，表示第 $i$ 个目标进行观测的时候采用的成像时间，如果不进行观测默认赋值为-1。
7. 实际观测开始时间： $ta_i$，表示第 $i$ 个目标进行观测的开始观测时间，如果不进行观测默认赋值为-1。
8. 观测顺序：$s_i$，表示第 $i$ 个目标被观测的顺序，如果不进行观测默认赋值为-1。

### 目标函数：
1. 最大化观测数目：f1=对于观测向量中非零的目标数进行求和操作即可得到；$f_1=-\sum_{i} x_i$

2. 最大化观测收益：f2=对于观测向量中非零的目标的观测收益（根据观测模式、成像模式）进行求和；$$f_2=-\sum_{i} x_i*price_{mode_i}$$
3. 最大化成像质量：f3=对于观测向量中非零的目标的图像质量（根据观测模式、成像模式）进行求和；$$f_2=-\sum_{i} x_i*quality_{mode_i}$$

### 约束条件：
#### 时间窗口约束
1. 最早观测开始时间不早于可见窗口开始时刻：$g1 = T_i.start - ta_i\le0$。
2. 最晚观测开始时间不晚于可见窗口结束时刻减去成像时长：$g2 = ta_i+dt_i-T_i.end \le0$
3. 下一次观测开始时间晚于当前观测结束时间：生成观测顺序序列的时候采用非重复序列生成的算法来进行顺序生成；

#### 成像区域相关约束：
1. 目标 $i$ 和目标 $j$ 的距离在成像宽度的一半以内：$g3 = dist(s_i,s_j)-\frac{amp_i}{2}\le0,\ if\ s_i+1==s_j$

### 输出需求
1. 获取每个目标的观测与否的状况；
2. 获取每个目标（特别是被观测的目标）的实际观测时间窗口开始时间以及采用的具体模式（是否侧视、采用的传感器模式）；

## 类似的问题

Minimize: TotalCommunicationTime

Subject to:

1. Each satellite is communicated exactly once:
$$
   \sum_{j \in GroundStations} x_{ij} = 1
$$
2. Each ground station communicates with exactly one satellite:
$$
   \sum_{i \in Satellites} x_{ij} = 1
$$
3. Time window constraints for satellite communication:
$$
   EarliestStart_i <= Arrival_j <= LatestEnd_i,  
$$
4. Total communication time constraint:
$$
   \sum_{i \in Satellites} communication_time_i * x_{ij} <= TotalCommunicationTime
$$
5. Binary decision variables:
$$
   x_{ij} = {0, 1}
$$