# refer:不同轨道类型LEO卫星轨道拟合及预报精度研究
import numpy as np
import pandas as pd

poly_data = pd.read_csv('/home/gwj/genetic/sat_data/track_simu_yaogan35.csv',header=None,index_col=None)
poly_data.columns = ['time_stamp','x','y','z','vx','vy','vz']

point_num = poly_data.shape[0]
x_data = np.linspace(0,10,point_num)
y_data = poly_data['x']

# 多项式阶数
poly_order = 6

# 拟合的数据点数
timepoints = 120

# 矩阵构造参数
n = poly_order
m = timepoints

# 拟合的数据点时间数据
x_arr = np.arange(timepoints)

# 数据点向量
Fx = y_data.iloc[:m].to_numpy()

coefficients = np.polynomial.chebyshev.chebfit(x_arr, Fx, poly_order)

# 代入预测点进行求解
pre_num = timepoints+600
pointx =np.arange(pre_num)
pred_y = []
pred_y = np.polynomial.chebyshev.chebval(pointx, coefficients)

# 对比原数据与预测的数据
print(pred_y)
import matplotlib.pyplot as plt
plt.figure(figsize=(16,8),dpi=200)
sunset_colors = ['#FFD28F', '#FFA64D', '#FF4D4D', '#FF33CC', '#B224EF']
plt.scatter(pointx,y_data.iloc[:pre_num],label='real data',color=sunset_colors[0],s=20)
plt.plot(pred_y,label='pre data', color=sunset_colors[1],linestyle='--')
plt.plot(y_data.iloc[:m],label='ori data',color = sunset_colors[2], linestyle='-')
plt.legend()
# plt.show()
plt.savefig("poly by chebyshev polynomial")

def calculate_rmse(predicted, actual):
    # 计算预测值与实际值之间的差的平方
    squared_errors = (predicted - actual) ** 2
    # 计算平均平方误差
    mean_squared_error = np.mean(squared_errors)
    # 计算 RMSE
    rmse = np.sqrt(mean_squared_error)
    return rmse

print(calculate_rmse(pred_y,y_data.iloc[:pre_num]))