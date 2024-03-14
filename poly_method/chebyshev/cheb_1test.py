# refer:不同轨道类型LEO卫星轨道拟合及预报精度研究
import numpy as np
import pandas as pd

poly_data = pd.read_csv('/home/gwj/genetic/sat_data/track_simu_yaogan35.csv',header=None,index_col=None)
poly_data.columns = ['time_stamp','x','y','z','vx','vy','vz']

point_num = poly_data.shape[0]
x_data = np.linspace(0,10,point_num)
y_data = poly_data['x']



# 因变量归一化
def normalize(t,th,tm):
    # 由于切比雪夫多项式自变量范围为，将对应区段因变量归一化到[-1,1]区间范围内
    return 2*(t-tm)/(th-tm)

# 切比雪夫多项式
def Tn(x,n):
    # T0=1
    # T1=x
    # T2=2x*T1-T0
    # ...
    # Tn=2x*T_{n-1}-T_{n-2}
    if n==0:
        return 1
    elif n==1:
        return x
    else:
        return 2*x*Tn(x,n-1)-Tn(x,n-2)

# 构造\tau对应的
def construct_T(x,n):
    T_list = [0 for i in range(n+1)]
    for i in range(n+1):
        T_list[i] = Tn(x,i)
    
    return T_list

# 多项式阶数
poly_order = 6

# 拟合的数据点数
timepoints = 32

# 矩阵构造参数
n = poly_order
m = timepoints

# 拟合的数据点时间数据
x_arr = np.arange(timepoints)

# 构造切比雪夫多项式矩阵
# m\times n的一个数组
T_matrix = np.array([[0 for i in range(n+1)] for j in range(m)])
for i in range(m):
    for j in range(n+1):
        T_matrix[i][j]= Tn(x_arr[i],j)

# 数据点向量
Fx = y_data.iloc[:m].to_numpy()

# 求解系数矩阵
# n+1\times m dot m\times n+1 = n+1\times n+1
TtT = np.dot(T_matrix.T,T_matrix)
# 转置之后矩阵仍为 n+1 x n+1 矩阵
TtTi = np.linalg.inv(TtT)
# n+1 \times m dot m = n+1\times 1
TtX = np.dot(T_matrix.T,Fx)
# n+1 x 1
Qx = np.dot(TtTi,TtX)
# print(Qx)

# 代入预测点进行求解
pointx =np.arange(timepoints)
pred_y = []
for pt in pointx:
    T_point = construct_T(pt,n)
    predx = Qx@T_point
    pred_y.append(predx)

# 对比原数据与预测的数据
print(pred_y)
import matplotlib.pyplot as plt
plt.plot(pred_y,label='pre data')
plt.plot(y_data.iloc[:m],label='ori data')
plt.legend()
plt.show()