import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Define the cosine function
def cosine_function(x, a, b, c):
    return a * np.cos(b * x + c)

# Generate some sample data
x_data = np.linspace(0, 10, 2000)  # Use limited data points
y_data = np.cos(x_data)*1803  # True cosine values

# Add some noise to the data (optional)
noise = np.random.normal(0, 0.1, len(y_data))
y_data += noise

# Fit the data to the cosine function
initial_guess = (1, 1, 1)  # Initial guess for parameters a, b, c
params, params_covariance = curve_fit(cosine_function, x_data, y_data, p0=initial_guess)

# Predict the curve using the fitted parameters
x_predict = np.linspace(0, 10, 100)  # Predict for a smoother curve
y_predict = cosine_function(x_predict, *params)

# Plot the original data and the predicted curve
plt.scatter(x_data, y_data, label='Data')
plt.plot(x_predict, y_predict, color='red', label='Predicted Curve')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Fitting Cosine Curve to Data')
plt.legend()
plt.show()
