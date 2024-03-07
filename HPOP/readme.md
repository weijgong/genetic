<!--
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-11-07 20:37:20
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-11-23 12:36:08
 * @FilePath: /gongweijing/High-Precision-Orbit-Propagator/readme.md
 * @Description: 
 * 
 * Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
-->

echo "Zx123456" | sudo -S git add .
echo "Zx123456" | sudo -S git commit -m "添加了HPOP算法"

# deploy
1. Download the gcc version == 4.8
```shell
sudo apt-get install aptitude
sudo aptitude install gcc-4.8
```
可能的报错：
- No candidate version found for gcc-4.8.Unable to apply some actions, aborting

解决方法：
<!-- 换源：不可行
```
打开source.list文件进行编辑：
sudo vim /etc/apt/sources.list
在末尾插入：
deb http://ppa.launchpad.net/ubuntu-toolchain-r/test/ubuntu xenial main 
deb-src http://ppa.launchpad.net/ubuntu-toolchain-r/test/ubuntu xenial main
激活：
sudo apt update
```
 -->
 <!-- wget https://ftp.gnu.org/gnu/gcc/gcc-4.8.0/gcc-4.8.0.tar.bz2 
version is not same with the tutorial
-->
<!-- wget http://10.88.162.17/software/gcc62/gmp-4.3.2.tar.bz2
这个源挂了 -->
<!-- # 安装mpc-0.8.1
wget http://10.88.162.17/software/gcc62/mpc-0.8.1.tar.gz
tar xvzf mpc-0.8.1.tar.gz 
没源
-->
<!-- # 安装mpc-1.0.1
wget https://mirrors.aliyun.com/gnu/mpc/mpc-1.0.1.tar.gz
tar xvzf mpc-1.0.1.tar.gz -->
首先，由于大部分apt源没有gcc482,所以，安装一系列的package：
```md
# 安装gcc-4.8.2
wget https://ftp.gnu.org/gnu/gcc/gcc-4.8.2/gcc-4.8.2.tar.bz2
tar xvjf gcc-4.8.2.tar.bz2
# 安装gmp-4.3.2
wget https://mirrors.aliyun.com/gnu/gmp/gmp-4.3.2.tar.bz2
tar xvjf gmp-4.3.2.tar.bz2 
# 安装mpc-0.8.1
ftp://gcc.gnu.org/pub/gcc/infrastructure/mpc-0.8.1.tar.gz
tar -zxvf mpc-0.8.1.tar.gz 

# mpfr-2.4.2
wget https://mirrors.aliyun.com/gnu/mpfr/mpfr-2.4.2.tar.bz2
tar xvjf mpfr-2.4.2.tar.bz2
# binutils-2.27
wget https://mirrors.aliyun.com/gnu/binutils/binutils-2.27.tar.bz2
tar xvjf binutils-2.27.tar.bz2
``` 
再次，对各个库进行编译：
```md
# gmp-432
cd gmp-4.3.2/
mkdir build && cd build
../configure --prefix=/opt/compiler/gcc-4.8.2
make -j20
sudo make install
# mpfr
cd /mpfr-2.4.2
mkdir build && cd build
../configure --prefix=/opt/compiler/gcc-4.8.2 --with-gmp=/opt/compiler/gcc-4.8.2
make -j20
sudo make install
# mpc-0.8.1
cd mpc-0.8.1/
mkdir build && cd build
../configure --prefix=/opt/compiler/gcc-4.8.2 --with-gmp=/opt/compiler/gcc-4.8.2 --with-mpfr=/opt/compiler/gcc-4.8.2
make -j20
sudo make install
# check out the lib
sudo make check

# gcc-4.8.2
export LD_LIBRARY_PATH=/opt/compiler/gcc-4.8.2/lib/:$LD_LIBRARY_PATH
cd gcc-4.8.2/
mkdir build && cd build
../configure --prefix=/opt/compiler/gcc-4.8.2 --with-gmp=/opt/compiler/gcc-4.8.2 --with-mpc=/opt/compiler/gcc-4.8.2 --with-mpfr=/opt/compiler/gcc-4.8.2 --enable-checking=release --enable-ld=yes --enable-gold=yes --disable-multilib
make -j20
make install
```
2. 

# source code
- SAT_VecMat：实现了矩阵的一系列运算的代码。例如，加减乘除、矩阵对角化、矩阵切片、取特定行元素、列元素等。
- SAT_Time：实现了Julian Date与Calender Date之间的日期转换以及日期输出的辅助函数。
- SAT_RefSys：实现了天体和地球参照系之间的转换,包含了多个坐标系之间的变换。
- SAT_Force：实现了基本的轨道动力学用到的参数计算方法。
- SAT_DE: 实现了常微分方程的数值积分方法，主要用于一些计算。
- SAT_Const：将一些常用的天体学常数进行存储及定义，例如地球半径、pi的数值等。

# data resource
- GGM03S.txt：地球重力场参数
- InitialState_120Sats.txt：初始的卫星轨道参数
- SatelliteStatesAll.txt：轨道初始状态、卫星参数、历元时间（Envisat）
- 