import pandas as pd
import numpy as np
from scipy.optimize import OptimizeResult
from sklearn.linear_model import LinearRegression
import scipy as sc
import keras as k
import tensorflow as tf
import matplotlib.pyplot as plt
import statsmodels.api as sm
from FE_Library import yearfrac
from scipy.stats import norm

# Load Time Series
quotes = pd.read_csv("DatasetPythonAss3.csv")
quotes["Date"] = pd.to_datetime(quotes["Date"], format="%d/%m/%Y")
quotes = quotes.set_index("Date")
AAPL = quotes["AAPL"]
SPX = quotes["SPX"]
# Plot Time Series
plt.figure(1)
plt.plot(AAPL, label='AAPL')
plt.title('AAPL Time-series')
plt.xlabel('Datas')
plt.ylabel('Prices')
plt.show()

plt.figure(2)
plt.plot(SPX, label='SPX')
plt.title('SPX Time-series')
plt.xlabel('Datas')
plt.ylabel('Prices')
plt.show()

# Compute Log_returns

vectorAAPL = quotes['AAPL'].to_numpy()
logReturnAAPL = np.divide(vectorAAPL[1:], vectorAAPL[:-1])
plt.figure(3)
plt.plot(logReturnAAPL, label='LogReturn AAPL')
plt.title('LogReturn AAPL Time-series')
plt.xlabel('Datas')
plt.ylabel('Prices')
plt.show()

vectorSPX = quotes['SPX'].to_numpy()
logReturnSPX = np.divide(vectorSPX[1:], vectorSPX[:-1])
plt.figure(4)
plt.plot(logReturnSPX, label='LogReturn SPX')
plt.title('LogReturn SPX Time-series')
plt.xlabel('Datas')
plt.ylabel('Prices')
plt.show()
# Regressions
model = LinearRegression()
logReturnSPX = logReturnSPX.reshape(-1, 1)
model.fit(logReturnSPX, logReturnAAPL)

plt.scatter(logReturnSPX, logReturnAAPL, label='Data')
plt.plot(logReturnSPX, model.predict(logReturnSPX), color='red')
plt.xlabel('X')
plt.ylabel('y')
plt.title('regression AAPL on SPX')
plt.legend()
plt.show()

slope = model.coef_[0]
print ("The slope is equal to", slope)
# YearFrac
first_element = quotes.index[0]
last_element = quotes.index[-1]

print("The first element is", first_element, "and the last one is", last_element)
fractionYear = yearfrac(first_element, last_element, 3)
print("The year frac between the first and the last element is",fractionYear)

# Interpolate
x = np.array([0, 1, 2, 3, 4]).reshape(-1, 1)
f = np.array([1, 2, 3.5, 4, 2])
model = LinearRegression()
model.fit(x, f)
point = np.array([[2.7]])
interp_val = model.predict(point)
print("The linear interpolation in", point, "is equal to",interp_val)
# Simulation of a standard normal random variable
np.random.seed(42)
standard_normal_samples = np.random.randn(100)
plt.hist(standard_normal_samples, bins=40, density=True, alpha=0.7, color='blue')
plt.title('Histogram of Standard Normal Distribution')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.show()

# show that the variance converges to 1 as the number of simulations increases
n = np.linspace(100, 10000, 100)
N = len(n)
variance = np.zeros(N)
for i in range(N):
    standard_normal_samples = np.random.randn(int(n[i]))
    variance[i] = np.var(standard_normal_samples)
plt.plot(n, variance)
plt.title('Variance of simulated standard gaussian')
plt.axhline(y=1, color='red', linestyle='--')
plt.xlabel('N')
plt.ylabel('Variance')
plt.show()

quantile_0_9_simulated = np.percentile(standard_normal_samples, 90)
quantile_0_9_real = norm.ppf(0.9)
CDF_at_quantile_0_9_real = norm.cdf(quantile_0_9_real)
CDF_at_quantile_0_9_simulated = norm.cdf(quantile_0_9_simulated)

print("The CDF evaluated with the real quantile 0.9 is", CDF_at_quantile_0_9_real)
print("The CDF evaluated with the simulated quantile 0.9 is", CDF_at_quantile_0_9_simulated)
print("Hence the absolute value of the difference between the two CDF is", abs(CDF_at_quantile_0_9_real-CDF_at_quantile_0_9_simulated))

# Minimization
def f(x):
    return (x[0]-3)**2+(x[1]-7)**2
x_y_min = sc.optimize.minimize(f,(0, 0))
print("The minimum of the function f is", x_y_min)
