import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# 定义余弦曲线用于进行数据点拟合
def cosine_function(x, a, b, c):
    return a * np.cos(b * x + c)

poly_data = pd.read_csv('/home/gwj/genetic/sat_data/track_simu_yaogan35.csv',header=None,index_col=None)
poly_data.columns = ['time_stamp','x','y','z','vx','vy','vz']
# 绘制各个曲线，确定拟合形式
# poly_data.plot()
# plt.savefig("/home/gwj/genetic/poly_method/origin_data.png")
# STK生成数据点的个数
point_num = poly_data.shape[0]
x_data = np.linspace(0,10,point_num)
y_data = poly_data['x']

degree = 5  # Degree of the polynomial
coefficients = np.polynomial.chebyshev.chebfit(x_data, y_data, degree)

x_pred = np.linspace(0,20,1000)
y_pred = np.polynomial.chebyshev.chebval(x_pred, coefficients)

# Plot the original data and the predicted curve
plt.scatter(x_data, y_data, label='Data')
plt.plot(x_pred, y_pred, color='red', label='Predicted Curve')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Fitting Polynomial using Chebyshev Basis Functions')
plt.legend()
plt.show()