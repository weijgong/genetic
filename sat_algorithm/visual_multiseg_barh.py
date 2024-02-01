import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import numpy as np

# 数据
segments = 5
categories = ['Category 1', 'Category 2', 'Category 3', 'Category 4', 'Category 5']
values = [[20, 15, 25, 10, 15],   # 第一个段
          [15, 20, 10, 25, 20],   # 第二个段
          [25, 10, 15, 20, 10],   # 第三个段
          [10, 25, 20, 15, 25],   # 第四个段
          [15, 15, 10, 30, 20]]   # 第五个段

# 绘图
fig, ax = plt.subplots()
y_pos = np.arange(len(categories))
left = np.zeros(len(categories))  # 初始化左边界为0

for i in range(segments):
    ax.barh(y_pos, values[i], align='center', left=left, label=f'Segment {i+1}')
    left += values[i]  # 更新左边界位置

ax.set_yticks(y_pos)
ax.set_yticklabels(categories)
ax.invert_yaxis()  # 反转 y 轴，使第一个类别显示在顶部
ax.set_xlabel('Values')
ax.set_title('Multi-segment Horizontal Bar Chart')
ax.legend()


plt.savefig("test_multisegment_barh.png")
