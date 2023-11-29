<!--
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-11-29 15:23:26
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-11-29 16:17:00
 * @FilePath: /gongweijing/nsga2/Project_Design.md
 * @Description: 记录项目主要框架
 * 
 * Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
-->
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

## STK simulate
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


# Question:
- what is J2000?
- 