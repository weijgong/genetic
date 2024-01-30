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

git add .
git commit -m "处理重载系统后配置及代码重写"
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
```
## 安装gnuplot适配初始绘图工具，后期完成总体py绘图配置可以不进行这步
```
sudo apt update
sudo apt install gnuplot
```
## 安装Python3.8
```
sudo apt install python3
apt install python3-pip
pip install matplotlib
```

# 文件内容
- time_trans.c ： 将UTC时间转换为世界秒(原子时)。
- LLA_Position.csv : 表示的是五个目标的lat,lon,alt.
- Sensor-to-target.csv : 表示的是五个目标的可见时间窗口。
- read_access.c ： 将可见时间窗口文件读入全局变量。
- 