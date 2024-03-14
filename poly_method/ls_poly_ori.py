import numpy as np
import matplotlib.pyplot as plt

# Generate some sample data
x_data = np.linspace(0, 10, 20)  # Use limited data points
y_data = np.cos(x_data)  # True cosine values

# Fit the data to a polynomial using least squares regression
degree = 5  # Degree of the polynomial
coefficients = np.polyfit(x_data, y_data, degree)

# Evaluate the polynomial at finer intervals for plotting
x_predict = np.linspace(0, 10, 100)
y_predict = np.polyval(coefficients, x_predict)

# Plot the original data and the fitted polynomial
plt.scatter(x_data, y_data, label='Data')
plt.plot(x_predict, y_predict, color='red', label='Fitted Polynomial')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Fitting Polynomial using Least Squares Regression')
plt.legend()
plt.show()
