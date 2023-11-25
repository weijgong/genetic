<!--
 * @Author: gongweijing 876887913@qq.com
 * @Date: 2023-11-24 11:02:30
 * @LastEditors: gongweijing 876887913@qq.com
 * @LastEditTime: 2023-11-24 15:10:52
 * @FilePath: /gongweijing/genetic/readme.md
 * @Description: Readme file for using test file
 * 
 * Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
-->
## run test-files
```md
# 运行程序
make
./tspmain --population 1000 --generation 10 --keep 100 --mutation 900 < ./test-files/test-30.txt
```

## debug test-files
```md
# debug 程序make
make
gdb tspmain
r --population 1000 --generation 10 --keep 100 --mutation 900 < ./test-files/test-30.txt

# 如果要观察 population 内容
break tsp-ga.cc:129
display population
```

## data 
```
30          the coordinate number of data
3 7 23      the coordinate X,Y,Z of one data
...
87 63 6
```