import numpy as np
import matplotlib.pyplot as plt

# Generate some sample data
x_data = np.linspace(0, 10, 20)  # Use limited data points
y_data = np.cos(x_data)  # True cosine values

# Fit the data to a polynomial using Chebyshev basis functions
degree = 5  # Degree of the polynomial
coefficients = np.polynomial.chebyshev.chebfit(x_data, y_data, degree)

# Evaluate the polynomial at finer intervals for plotting
x_predict = np.linspace(0, 10, 100)
y_predict = np.polynomial.chebyshev.chebval(x_predict, coefficients)

# Plot the original data and the fitted polynomial
plt.scatter(x_data, y_data, label='Data')
plt.plot(x_predict, y_predict, color='red', label='Fitted Polynomial')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Fitting Polynomial using Chebyshev Basis Functions')
plt.legend()
plt.show()
