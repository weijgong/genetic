<!--
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-12-02 00:52:44
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-12-03 00:15:23
 * @FilePath: /root/genetic/sat_algorithm/readme.md
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
-->

为了方便下次进行重装系统配置环境，特此记录配置过程：
## 配置git
```
sudo su
git config --global user.email "876887913@qq.com"
git config --global user.name "weijgong"

sudo su
echo "Zx123456" | sudo -S git add .
echo "Zx123456" | sudo -S git commit -m "完成可视化流图代码，部分重规划代码"


echo "Zx123456" | sudo -S git add .
echo "Zx123456" | sudo -S git commit -m "修复了tournament中的浅拷贝问题，考虑未来采用vector替代E*"

echo "Zx123456" | sudo -S git add .
echo "Zx123456" | sudo -S git commit -m "暂时进行代码完全重构，将原本的指针数组进行修改，防止因指针产生的复杂问题。"


echo "Zx123456" | sudo -S git add .
echo "Zx123456" | sudo -S git commit -m "进行了重构"

echo "Zx123456" | sudo -S git add .
echo "Zx123456" | sudo -S git commit -m "修复了mktime过程中的空指针的bug"
git push origin nsga-II
```

## sudo更新
```
sudo apt update
sudo apt upgrade
```

## 安装gcc、g++等配件
```
sudo apt update
sudo apt install g++
sudo apt install build-essential
sudo apt install gdb
sudo apt install valgrind
```
## 安装gnuplot适配初始绘图工具，后期完成总体py绘图配置可以不进行这步
```
sudo apt update
sudo apt install gnuplot
```
## 安装Python3.8
```
sudo -S sudo apt install python3
sudo -S sudo apt install python3-pip
sudo -S pip3 install matplotlib
sudo -S sudo apt-get install python3-tk
```

# 文件内容
- time_trans.c ： 将UTC时间转换为世界秒(原子时)。
- LLA_Position.csv : 表示的是五个目标的lat,lon,alt.
- Sensor-to-target.csv : 表示的是五个目标的可见时间窗口。
- read_access.c ： 将可见时间窗口文件读入全局变量。
- 