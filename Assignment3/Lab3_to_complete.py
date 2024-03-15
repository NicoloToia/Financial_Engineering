import pandas as pd
import numpy as np
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
logReturnAAPL = np.log(1 + AAPL.pct_change())
logReturnSPX = np.log(1 + SPX.pct_change())

# delete the first element (NaN)
logReturnAAPL = logReturnAAPL[1:]
logReturnSPX = logReturnSPX[1:]

# plot log returns

plt.figure(3)
plt.plot(logReturnAAPL, label='LogReturn AAPL')
plt.title('LogReturn AAPL Time-series')
plt.xlabel('Datas')
plt.ylabel('Prices')
plt.show()

plt.figure(4)
plt.plot(logReturnSPX, label='LogReturn SPX')
plt.title('LogReturn SPX Time-series')
plt.xlabel('Datas')
plt.ylabel('Prices')
plt.show()

# Regressions
model = LinearRegression()
model.fit(logReturnSPX.values.reshape(-1, 1), logReturnAAPL.values.reshape(-1, 1))

plt.scatter(logReturnSPX, logReturnAAPL, label='Data')
plt.plot(logReturnSPX, model.predict(logReturnSPX.values.reshape(-1, 1)), label='Linear regression', color='red')
plt.xlabel('X')
plt.ylabel('y')
plt.legend()
plt.show()

slope = model.coef_[0]
print (slope)
# YearFrac
first_element = quotes.index[0]
last_element = quotes.index[-1]
print(first_element, last_element)
fractionYear = yearfrac(first_element, last_element, 3)
print(fractionYear)

# Interpolate
x = np.array([0, 1, 2, 3, 4]).reshape(-1, 1)
f = np.array([1, 2, 3.5, 4, 2])
model = LinearRegression()
model.fit(x, f)
point = np.array([[2.7]])
interp_val = model.predict(point)
print (interp_val)
# Simulation
np.random.seed(42)
#n = np.linspace(100, 1000, 100)
#N = len(n)
#standard_normal_samples = np.array([])
#variance = np.array([])
#for i in range(0, N-1):
    #standard_normal_samples[i] = np.random.randn((i+1)*100)
    #variance[i] = np.var(standard_normal_samples[i])
standard_normal_samples= np.random.randn(100)
variance100= np.var(standard_normal_samples)
plt.hist(standard_normal_samples, bins=40, density=True, alpha=0.7, color='blue')
plt.title('Histogram of Standard Normal Distribution')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.show()

standard_normal_samples= np.random.randn(300)
variance300= np.var(standard_normal_samples)
plt.hist(standard_normal_samples, bins=40, density=True, alpha=0.7, color='blue')
plt.title('Histogram of Standard Normal Distribution')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.show()

standard_normal_samples= np.random.randn(500)
variance500= np.var(standard_normal_samples)
plt.hist(standard_normal_samples, bins=40, density=True, alpha=0.7, color='blue')
plt.title('Histogram of Standard Normal Distribution')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.show()

standard_normal_samples= np.random.randn(800)
variance800= np.var(standard_normal_samples)
plt.hist(standard_normal_samples, bins=40, density=True, alpha=0.7, color='blue')
plt.title('Histogram of Standard Normal Distribution')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.show()

standard_normal_samples= np.random.randn(1000)
variance1000= np.var(standard_normal_samples)
plt.hist(standard_normal_samples, bins=40, density=True, alpha=0.7, color='blue')
plt.title('Histogram of Standard Normal Distribution')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.show()

n = ([100, 300, 500, 800, 1000])
variance = ([variance100, variance300,variance500, variance800, variance1000])
plt.plot(n, variance)
plt.show()

quantile_0_9 = np.percentile(standard_normal_samples, 90)
CDF_at_quantile_0_9 = norm.cdf(quantile_0_9)

# Minimization


