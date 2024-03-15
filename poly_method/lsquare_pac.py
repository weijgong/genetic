# refer:不同轨道类型LEO卫星轨道拟合及预报精度研究
import numpy as np
import pandas as pd

poly_data = pd.read_csv('/home/gwj/genetic/sat_data/track_simu_yaogan35.csv',header=None,index_col=None)
poly_data.columns = ['time_stamp','x','y','z','vx','vy','vz']

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
Fx  = poly_data['x' ].iloc[:m].to_numpy()
Fy  = poly_data['y' ].iloc[:m].to_numpy()
Fz  = poly_data['z' ].iloc[:m].to_numpy()
Fvx = poly_data['vx'].iloc[:m].to_numpy()
Fvy = poly_data['vy'].iloc[:m].to_numpy()
Fvz = poly_data['vz'].iloc[:m].to_numpy()
import time
st = time.time()
# 拟合系数矩阵
coef_x  = np.polyfit(x_arr, Fx , poly_order)
coef_y  = np.polyfit(x_arr, Fy , poly_order)
coef_z  = np.polyfit(x_arr, Fz , poly_order)
coef_vx = np.polyfit(x_arr, Fvx, poly_order)
coef_vy = np.polyfit(x_arr, Fvy, poly_order)
coef_vz = np.polyfit(x_arr, Fvz, poly_order)

# 代入预测点进行求解
pre_num = timepoints+600
pointx =np.arange(pre_num)
pred_x  = np.polyval(coef_x , pointx)
pred_y  = np.polyval(coef_y , pointx)
pred_z  = np.polyval(coef_z , pointx)
pred_vx = np.polyval(coef_vx, pointx)
pred_vy = np.polyval(coef_vy, pointx)
pred_vz = np.polyval(coef_vz, pointx)
print(time.time()-st)
# 对比原数据与预测的数据
# print(pred_y)
import matplotlib.pyplot as plt
plt.figure(figsize=(16,8),dpi=200)
sunset_colors = ['#FFD28F', '#FFA64D', '#FF4D4D', '#FF33CC', '#B224EF','#fd8a27']
oceanic_colors = ['#8ED6FF', '#5D9CEC', '#265A84', '#48CFAD', '#16A085','#a5d6f8']
plt.scatter(pointx,poly_data['x' ].iloc[:pre_num],label='real data x ',color=oceanic_colors[0],s=20)
plt.scatter(pointx,poly_data['y' ].iloc[:pre_num],label='real data y ',color=oceanic_colors[1],s=20)
plt.scatter(pointx,poly_data['z' ].iloc[:pre_num],label='real data z ',color=oceanic_colors[2],s=20)


plt.plot(pred_x ,label='pre data x', color=sunset_colors[0],linestyle='--')
plt.plot(pred_y ,label='pre data y', color=sunset_colors[1],linestyle='--')
plt.plot(pred_z ,label='pre data z', color=sunset_colors[2],linestyle='--')

plt.legend()
plt.title('polynomial of degree n predict xyz')
plt.xlabel('time points num')
plt.ylabel('km')
plt.grid(alpha=0.4,linestyle='--')
# plt.show()
plt.savefig("poly by polynomial of degree n xyz")

'''
plt.figure(figsize=(16,8),dpi=200)
plt.scatter(pointx,poly_data['vx'].iloc[:pre_num],label='real data vx',color=sunset_colors[3],s=20)
plt.scatter(pointx,poly_data['vy'].iloc[:pre_num],label='real data vy',color=sunset_colors[4],s=20)
plt.scatter(pointx,poly_data['vz'].iloc[:pre_num],label='real data vz',color=sunset_colors[5],s=20)

plt.plot(pred_vx,label='pre datavx', color=oceanic_colors[3],linestyle='--')
plt.plot(pred_vy,label='pre datavy', color=oceanic_colors[4],linestyle='--')
plt.plot(pred_vz,label='pre datavz', color=oceanic_colors[5],linestyle='--')
plt.legend()
# plt.show()
plt.savefig("poly by polynomial vx vy vz")
'''


def calculate_rmse(predicted, actual):
    mean_squared_error = 0
    for i in range(len(predicted)):
        squared_errors = (predicted[i] - actual[i]) ** 2
        mean_squared_error += np.mean(squared_errors)
    rmse = np.sqrt(mean_squared_error)
    return rmse

predicted = [pred_x,pred_y,pred_z]
y_data = [poly_data['x' ].iloc[:pre_num],poly_data['y' ].iloc[:pre_num],poly_data['z' ].iloc[:pre_num]]
print(calculate_rmse(predicted,y_data))
# 1.237192