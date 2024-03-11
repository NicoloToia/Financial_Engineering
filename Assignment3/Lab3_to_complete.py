import pandas as pd
import numpy as np
import scipy as sc
import keras as k
import tensorflow as tf
import matplotlib.pyplot as plt
import statsmodels.api as sm
from FE_Library import yearfrac

# Load Time Series
quotes = pd.read_csv("DatasetPythonAss3.csv")
quotes["Date"] = pd.to_datetime(quotes["Date"], format="%d/%m/%Y")
quotes = quotes.set_index("Date")
AAPL = quotes["AAPL"]
SPX = quotes["SPX"]
# Plot Time Series


# Compute Log_returns

# Regressions

# YearFrac

# Interpolate

# Simulation

# Minimization


