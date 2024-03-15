# import the necessary libraries
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

## Exercise 1

# Load Time Series
quotes = pd.read_csv("DatasetPythonAss3.csv")
quotes["Date"] = pd.to_datetime(quotes["Date"], format="%d/%m/%Y")
quotes = quotes.set_index("Date")
AAPL = quotes["AAPL"]
SPX = quotes["SPX"]

## Exercise 2

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

## Exercise 3

# Regression
model = LinearRegression()
model.fit(logReturnSPX.values.reshape(-1, 1), logReturnAAPL.values.reshape(-1, 1))

plt.scatter(logReturnSPX, logReturnAAPL, label='Data')
plt.plot(logReturnSPX, model.predict(logReturnSPX.values.reshape(-1, 1)), label='Linear regression', color='red')
plt.xlabel('X')
plt.ylabel('y')
plt.title('regression AAPL on SPX')
plt.legend()
plt.show()

slope = model.coef_[0]
print (f"The slope of the linear regression is {slope[0]:.4f}")

## Exercise 4

# YearFrac
first_element = quotes.index[0]
last_element = quotes.index[-1]

fractionYear = yearfrac(first_element, last_element, 3)
print(f"The fraction of year between the first and the last element is {fractionYear:.4f}")

## Exercise 5

# Interpolate
fx = [1, 2, 3.5, 4, 2]
x = [0, 1, 2, 3, 4]

# linearn interpolation
x_interp = 2.7
y_interp = np.interp(x_interp, x, fx)

print(f"The linear interpolation of {x_interp} is {y_interp:.4f}")

## Exercise 6

# Simulation of a standard normal random variable
np.random.seed(42)

# show that the variance converges to 1 as the number of simulations increases
n = np.linspace(100, 10000, 100)
variances = [
    np.var(np.random.randn(int(i)))
    for i in n
]

plt.plot(n, variances)
plt.title('Variance of simulated standard gaussian')
plt.axhline(y=1, color='red', linestyle='--')
plt.xlabel('N')
plt.ylabel('Variance')
plt.show()

# show that the gaussian CDF of the simulated data gives the correct result
z_sample = np.random.randn(10000)
quantile_90 = np.percentile(z_sample, 90)
CDF_sample = norm.cdf(quantile_90)

# print the result
print(f"The CDF evaluated with the simulated quantile 0.9 is {CDF_sample:.4f}")

## Exercise 7

# minimize the function (x-3)^2 + (y-7)^2
minimization = sc.optimize.minimize(lambda x: (x[0]-3)**2 + (x[1]-7)**2, [0, 0])

# print the coordinates
print(f"The minimum is found in ({minimization.x[0]:.4f}, {minimization.x[1]:.4f})")

