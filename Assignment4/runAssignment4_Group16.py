# import needed modules
from contextlib import redirect_stdout
import pandas as pd
import numpy as np
from scipy.stats import t
from scipy.stats import norm
from datetime import datetime
from statsmodels.multivariate.pca import PCA
# import custom functions and modules
from Assignment4_lib import HSMeasurements, bootstrapStatistical, WHSMeasurements, \
    PrincCompAnalysis, plausibilityCheck

# overwrite the print function to write to a file
def print(*args, **kwargs):
    # execute the normal print function
    __builtins__.print(*args, **kwargs)
    # write the output to a file
    with open('output.txt', 'a') as f:
        with redirect_stdout(f):
            __builtins__.print(*args, **kwargs)

# flush the file
with open('output.txt', 'w') as f:
    pass

# fix the random seed
np.random.seed(42)

# load the indexes as dictionary of {ticker: name}
with open('data/_indexes.csv', 'r') as f:
    # skip the first line
    indexes = {
        line.split(',')[1]: line.split(',')[2].strip()
        for line in f.readlines()[1:]
    }

# load the returns dataset as a dataframe
EuroStoxx50 = pd.read_csv('data/EUROSTOXX50_Dataset.csv', sep=',')
EuroStoxx50 = EuroStoxx50.set_index(pd.DatetimeIndex(EuroStoxx50['Date']))
EuroStoxx50 = EuroStoxx50.drop('Date', axis=1)

# get the head and describe the dataset
EuroStoxx50.head()
EuroStoxx50.describe()

# Data Cleaning and Preprocessing

# delete columns of only missing values and fill others with the previous value
EuroStoxx50 = EuroStoxx50.dropna(axis=1, how='all')
EuroStoxx50 = EuroStoxx50.ffill()

# create the log-returns dataframe
returns = np.log(EuroStoxx50/EuroStoxx50.shift(1))
returns = returns.dropna(axis=0, how='all')

returns.head()

## POINT 0: Variance-Covariance Method for VaR and ES in a linear portfolio
# equally weighted portfolio made up of
# - Adidas (ADSGn.DE)
# - Allianz (ALVG.DE)
# - Munich Re (MUVGn.DE)
# - L'Oreal (OREP.PA)

# set the parameters
nu_0 = 4
alpha_0 = 0.99
notional_0 = 15 * 10**6
# weights of the portfolio
w_0 = np.ones(4) / 4

# select the returns of the portfolio
df_0 = returns[['ADSGn.DE', 'ALVG.DE', 'MUVGn.DE', 'OREP.PA']]
valuation_date_0 = datetime(2020, 2, 20)
# only use data prior to the valuation date and no older than 5 years
df_0 = df_0[df_0.index <= valuation_date_0]
df_0 = df_0[df_0.index > valuation_date_0 - pd.DateOffset(years=5)]

# compute mean vector and covariance matrix
mu_0 = df_0.mean()
Sigma_0 = df_0.cov()

# Daily VaR

t_alpha = t.ppf(alpha_0, nu_0)

VaR_0 = mu_0 @ w_0 + np.sqrt(w_0 @ Sigma_0 @ w_0) * t_alpha

print(f"""
 >--- Point 0: t-student VaR ---<
 The daily VaR at 99% confidence level is: {VaR_0:.2%}
 The daily VaR at 99% confidence level is: {VaR_0 * notional_0:.2f} EUR
""")

# Daily ES

ES_std = (nu_0 + t_alpha**2) / (nu_0 - 1) * (t.pdf(t_alpha, nu_0) / (1 - alpha_0))

ES_0 = mu_0 @ w_0 + np.sqrt(w_0 @ Sigma_0 @ w_0) * ES_std

print(f"""
 >--- Point 0: t-student ES ---<
 The daily ES at 99% confidence level is: {ES_0:.2%}
 The daily ES at 99% confidence level is: {ES_0 * notional_0:.2f} EUR
""")

## Point 1: Historical Simulation, Statistical Bootstrap and PCA approach for VaR and ES in a
# linear portfolio

# set the parameters and data for the exercise
alpha_1 = 0.95
lmd_1 = 0.94 # lambda is a reserved keyword
valuation_date_1 = datetime(2019, 3, 20)

# select only the relevant returns up to the valuation date and no older than 5 years
df_1 = returns[returns.index <= valuation_date_1]
df_1 = df_1[df_1.index > valuation_date_1 - pd.DateOffset(years=5)]

## Point 1.1: Historical Simulation
# Portfolio made up of
# = 25K shares of Total (TTEF.PA)
# = 20K shares of AXA (AXAF.PA)
# - 20K shares of Sanofi (SASY.PA)
# - 10K shares of Volkswagen (VOWG_p.DE)

# select the relevant indexes
df_1_1 = df_1[['TTEF.PA', 'AXAF.PA', 'SASY.PA', 'VOWG_p.DE']]
# compute the value at valuation date
val_Total = 25_000 * EuroStoxx50.loc[valuation_date_1]['TTEF.PA']
val_AXA = 20_000 * EuroStoxx50.loc[valuation_date_1]['AXAF.PA']
val_Sanofi = 20_000 * EuroStoxx50.loc[valuation_date_1]['SASY.PA']
val_VW = 10_000 * EuroStoxx50.loc[valuation_date_1]['VOWG_p.DE']
# reconstruct the value and compute the weights
V_t_1_1 = val_Total + val_AXA + val_Sanofi + val_VW
w_1_1 = np.array([val_Total, val_AXA, val_Sanofi, val_VW]) / V_t_1_1

# compute the historical VaR and ES with HSMeasurements
ES_HS, VaR_HS = HSMeasurements(df_1_1, alpha_1, w_1_1, V_t_1_1, 1)

print(f"""
 >--- Historical Simulation ---<
 alpha = {alpha_1} V_t = {V_t_1_1:.2f} EUR
    -> Daily VaR: {VaR_HS:>10.2f} EUR
    -> Daily ES: {ES_HS:>11.2f} EUR
""")

## Point 1.1: Statistical Bootstrap

# sample with bootstrapStatistical
sample_1_1 = bootstrapStatistical(200, df_1_1)
# compute the same as with HSMeasurements
ES_BS, VaR_BS = HSMeasurements(sample_1_1, alpha_1, w_1_1, V_t_1_1, 1)

print(f"""
 >--- Bootstrap ---<
 alpha = {alpha_1} V_t = {V_t_1_1:.2f} EUR
    -> Daily VaR: {VaR_BS:>10.2f} EUR
    -> Daily ES: {ES_BS:>11.2f} EUR
""")

# Plausibility Check for 1.1

# compute the plausibility check
VaR_PC_1_1 = plausibilityCheck(df_1_1, w_1_1, alpha_1, V_t_1_1, 1)

print(f"""
 >--- Plausibility Check ---<
    -> Thumb Rule: {VaR_PC_1_1:12.2f} EUR
    -> Against the HS: {VaR_HS:.2f} EUR
    -> Against the BS: {VaR_BS:.2f} EUR
""")

## Point 1.2: Weighted Historical Simulation
# equally weighted portfolio made up of
# - Adidas (ADSGn.DE)
# - Airbus (AIR.PA)
# - BBVA (BBVA.MC)
# - BMW (BMWG.DE)
# - Deutsche Telekom (DTEGn.DE)

# select the relevant indexes
df_1_2 = df_1[['ADSGn.DE', 'AIR.PA', 'BBVA.MC', 'BMWG.DE', 'DTEGn.DE']]
V_t_1_2 = 1
w_1_2 = np.ones(5) / 5

# compute the historical VaR and ES with WHSMeasurements
ES_WHS, VaR_WHS = WHSMeasurements(df_1_2, alpha_1, lmd_1, w_1_2, V_t_1_2, 1)

print(f"""
 >--- Weighted Historical Simulation ---<
 alpha = {alpha_1} lambda = {lmd_1}
    -> Daily VaR: {VaR_WHS:>10.2%}
    -> Daily ES: {ES_WHS:>11.2%}
""")

# Plausibility Check for 1.2
VaR_PC_1_2 = plausibilityCheck(df_1_2, w_1_2, alpha_1, V_t_1_2, 1)

print(f"""
 >--- Plausibility Check ---<
    -> Thumb Rule: {VaR_PC_1_2:>6.2%}
    -> Against the WHS: {VaR_WHS:>.2%}
""")

## Point 1.3: Principal Component Analysis
# We take the first 18 returns of the dataset and apply PCA for n = 1,...,5

# select the relevant indexes
df_1_3 = df_1.dropna(axis=1, how='any').iloc[:, :18]
Vt_1_3 = 1
w_1_3 = np.ones(18) / 18

# compute the yearly mean and covariance matrix
mu_1_3 = 256 * df_1_3.mean()
Sigma_1_3 = 256 * df_1_3.cov()

for i in range(1,6):
    # compute the VaR and ES with PrincCompAnalysis
    ES_PCA, VaR_PCA = PrincCompAnalysis(Sigma_1_3, mu_1_3, w_1_3, 10/256, alpha_1, i, Vt_1_3)

    print(f"""
 >--- Principal Component Analysis (n={i}) ---<
    -> Daily VaR: {VaR_PCA:>10.2%}
    -> Daily ES: {ES_PCA:>11.2%}
    """)

# Plausibility Check for 1.3
VaR_PC_1_3 = plausibilityCheck(df_1_3, w_1_3, alpha_1, Vt_1_3, 10)

print(f"""
 >--- Plausibility Check ---<
    -> Thumb Rule: {VaR_PC_1_3:>6.2%}
    -> Against the PCA: {VaR_PCA:>.2%}
""")
