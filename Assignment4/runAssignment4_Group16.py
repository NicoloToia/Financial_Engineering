# import needed modules
from contextlib import redirect_stdout
import pandas as pd
import numpy as np
from datetime import datetime
from scipy.stats import norm
# import custom functions and modules
from Assignment4_lib import AnalyticalNormalMeasures, HSMeasurements, bootstrapStatistical, \
    WHSMeasurements, PrincCompAnalysis, plausibilityCheck, FullMonteCarloVaR, DeltaNormalVaR, \
    DeltaGammaVaR, blackScholesCall

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

ES_0, VaR_0 = AnalyticalNormalMeasures(alpha_0, nu_0, w_0, notional_0, 1, df_0)

VaR_PC_0 = plausibilityCheck(df_0, w_0, alpha_0, notional_0, 1)

print("""
        >--- POINT 0 ---<
    """)

print(f"""
 >--- Point 0: t-student VaR ---<
 The daily VaR at 99% confidence level is: {VaR_0:.2f}
 The daily ES at 99% confidence level is: {ES_0:.2f}
 The plausibility check gives: {VaR_PC_0:.2f}
""")

## Point 1: Historical Simulation, Statistical Bootstrap and PCA approach for VaR and ES in a
# linear portfolio

# set the parameters and data for the exercise
alpha_1 = 0.95
lmd_1 = 0.95 # lambda is a reserved keyword
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
print("""
        >--- POINT 1 ---<
    """)
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
    -> Daily VaR: {VaR_WHS:>10.6%}
    -> Daily ES: {ES_WHS:>11.6%}
""")

# Plausibility Check for 1.2
VaR_PC_1_2 = plausibilityCheck(df_1_2, w_1_2, alpha_1, V_t_1_2, 1)

print(f"""
 >--- Plausibility Check ---<
    -> Thumb Rule: {VaR_PC_1_2:>6.6%}
    -> Against the WHS: {VaR_WHS:>.6%}
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
    ES_PCA, VaR_PCA, exp_var = PrincCompAnalysis(Sigma_1_3, mu_1_3, w_1_3, 10/256, alpha_1, i, Vt_1_3)

    print(f"""
 >--- Principal Component Analysis (n={i}, {exp_var:.2%}) ---<
    -> Daily VaR: {VaR_PCA:>10.6%}
    -> Daily ES: {ES_PCA:>11.6%}
    """)

# Plausibility Check for 1.3
VaR_PC_1_3 = plausibilityCheck(df_1_3, w_1_3, alpha_1, Vt_1_3, 10)

print(f"""
 >--- Plausibility Check ---<
    -> Thumb Rule: {VaR_PC_1_3:>.6%}
    -> Against the PCA: {VaR_PCA:>.6%}
""")

## Point 2: Full-Monte Carlo and Delta Normal Approach for VaR
print("""
        >--- POINT 2 ---<
    """)

valuation_date_2 = datetime(2017, 1, 16)
num_days_in_year = 256

# retrieve the relevant returns (up to 2 years before the valuation date)
df_2 = returns[returns.index <= valuation_date_2]
df_2 = df_2[df_2.index > valuation_date_2 - pd.DateOffset(years=2)]
df_2 = df_2['BMWG.DE']

# call parameters (yearly data)
expiry_date_2 = datetime(2017, 4, 18)
K = 25 # strike price
sigma = 15.4 / 100 # volatility
delta = 3.1 / 100 # dividend yield
r = 0.5 / 100 # risk-free rate

ttm = (expiry_date_2 - valuation_date_2).days / num_days_in_year

# set the parameters and data for the exercise
alpha_2 = 0.95
lmd_2 = 0.95
H = 10 / num_days_in_year

# deduce the number of shares
total_val_BMW = 1_186_680
S_t = EuroStoxx50.loc[valuation_date_2]['BMWG.DE']
num_shares_BMW = total_val_BMW / S_t
num_calls = num_shares_BMW

# compute the VaR with Monte Carlo
VaR_MC = FullMonteCarloVaR(df_2, num_shares_BMW, num_calls, S_t, K, r, delta, sigma, ttm, H, alpha_2, lmd_2, num_days_in_year)

print(f"""
 >--- Full Monte Carlo VaR ---<
    -> 10-days VaR: {VaR_MC:.2f}
""")

VaR_DN = DeltaNormalVaR(df_2, num_shares_BMW, num_calls, S_t, K, r, delta, sigma, ttm, H, alpha_2, lmd_2, num_days_in_year)

print(f""" 
 >--- Delta Normal VaR ---<
    -> 10-days VaR: {VaR_DN:.6f}
""")

VaR_DG = DeltaGammaVaR(df_2, num_shares_BMW, num_calls, S_t, K, r, delta, sigma, ttm, H, alpha_2, lmd_2, num_days_in_year)

print(f"""
 >--- Delta Gamma VaR ---<
    -> 10-days VaR: {VaR_DG:.6f}
""")


# Point 3: Clicquet option
print("""
        >--- POINT 3 ---<
    """)

# problem parameters
valuation_date_3 = datetime(2008, 2, 19)
notional_3 = 30 * 10**6
L = 0.99
sigma_3 = 20 / 100

# load the discounts and probabilities for the relevant dates
P_ISP = pd.read_csv('data/P_ISP.csv', sep=',', index_col=0, parse_dates=True)
DF = pd.read_csv('data/discountsCDS.csv', sep=',', index_col=0, parse_dates=True)

# MC simulation to compute the price of the Clicquet option
N_sim = 10**6
N_steps = len(DF)

# first MC to simulate the time of default tau

U = np.random.uniform(size=(N_sim, ))

# inverse survival function
default_idx = [
    len(P_ISP[P_ISP['survProb'] > u])
    for u in U
]

g = np.random.normal(size=(N_sim, N_steps))

S_t = np.zeros((N_sim, N_steps+1))

# set the initial value
S_t[:, 0] = 1

for i in range(N_steps):

    # compute the forward discount factor
    fwd_DF = DF['discount'].iloc[i] / DF['discount'].iloc[i-1] if i > 0 else DF['discount'].iloc[0]

    # compute the drift, diffusion and time step
    dt = (DF.index[i] - DF.index[i-1]).days / 365 if i > 0 else (DF.index[i] - valuation_date_3).days / 365

    drift = - sigma_3**2 / 2 * dt
    diffusion = sigma_3 * np.sqrt(dt)

    S_t[:, i+1] = S_t[:, i] * 1/fwd_DF * np.exp(drift + diffusion * g[:, i])

# plt.figure()
# for i in range(N_sim):
#     plt.plot(S_t[i, :])
# plt.show()

# compute the coupon payments as max(L * S_t - S_{t-1}, 0)

coupons = np.maximum(L * S_t[:, 1:] - S_t[:, :-1], 0)

R = 0.4

# loop over the paths
NPV = np.zeros(N_sim)
for i in range(N_sim):
    # no default
    if default_idx[i] == N_steps:
        NPV[i] = (coupons[i,:] * DF['discount'].values).sum()
    else:
        NPV[i] = (coupons[i, :default_idx[i]] * DF['discount'].values[:default_idx[i]]).sum() + \
            R * (coupons[i, default_idx[i]:] * DF['discount'].values[default_idx[i]:]).sum()

# compute the variance for the IC
std_dev = NPV.std()

NPV = NPV.mean()

alpha = 0.95
z_99 = norm.ppf(alpha)

# compute the IC with the notional and the confidence level
IC = [NPV - z_99 * std_dev / np.sqrt(N_sim), NPV + z_99 * std_dev / np.sqrt(N_sim)]

IC[0], IC[1] = IC[0] * notional_3, IC[1] * notional_3

print(f"""
 >--- Clicquet Option (MC) ---<
    The price of the Clicquet option is: {NPV:.8f} EUR
    Notional: {NPV*notional_3 :.8f}EUR
    The {alpha:.2%} confidence interval is: {IC[0]:.8f} - {IC[1]:.8f} EUR
""")

## closed formula for the Clicquet option

s = 0
S_0 = 1

for i in range(N_steps):
    
    P = P_ISP['survProb'].iloc[i]
    prevProb = P_ISP['survProb'].iloc[i-1] if i > 0 else 1

    # compute the forward discount factor
    fwd_DF = DF['discount'].iloc[i] / DF['discount'].iloc[i-1] if i > 0 else DF['discount'].iloc[0]

    # compute the yearfrac
    yf = (DF.index[i] - DF.index[i-1]).days / 365 if i > 0 else (DF.index[i] - valuation_date_3).days / 365

    C_t = blackScholesCall(L / fwd_DF, 1, 0, 0, sigma_3, yf)

    first_term = fwd_DF * S_0 * P * C_t

    second_term = 0

    for j in range(i, N_steps):

        fwd = DF['discount'].iloc[j] / DF['discount'].iloc[j-1] if j > 0 else DF['discount'].iloc[j]
        yf = (DF.index[j] - DF.index[j-1]).days / 365 if j > 0 else (DF.index[j] - valuation_date_3).days / 365

        C_t = blackScholesCall(L / fwd, 1, 0, 0, sigma_3, yf)

        second_term += C_t

    second_term = R * S_0 * (prevProb - P) * second_term

    s += first_term + second_term

print(f"""
 >--- Closed Formula Clicquet Option ---<
    The price of the Clicquet option is: {s:.8f} EUR
    Notional: {s*notional_3:.8f} EUR
""")

# The total price considering the notional
print(f"""
    >--- Closed Formula Clicquet Option price for 30mln notional ---<
      The total price of the Clicquet option is: {s * notional_3:.2f} EUR)
""")