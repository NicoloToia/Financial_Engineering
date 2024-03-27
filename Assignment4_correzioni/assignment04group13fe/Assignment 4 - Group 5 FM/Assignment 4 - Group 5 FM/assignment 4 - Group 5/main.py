from utilities import *

# Importing the data
data = pd.read_csv("EUROSTOXX50_Dataset.csv", index_col=0)
index = ["ADSGn.DE", "ALVG.DE", "MUVGn.DE", "OREP.PA"]
data1 = data[index].copy()

# Adding previous day value in case of missing price
data1 = data1.fillna(method="ffill")

# Indexes as datetime type
data1.index = pd.to_datetime(data1.index)

# Taking the return in the last 5 years
current_date = datetime.datetime(2020, 2, 20)
d = 365 * 5 + 1  # There is one leap year
delta = datetime.timedelta(days=d)
past = current_date - delta
data1 = data1.loc[past:current_date]

# Equally weighted portfolio
w = np.array([0.25, 0.25, 0.25, 0.25])

# Computing log-returns
# shift(0) is the "current" return, while shift(1) is the previous one
log_returns = np.log(data1.shift(0) / data1.shift(1))

# Eliminating NaN from the dataframe
log_returns = log_returns.dropna()

# Setting parameters
alpha = 0.99
notional = 15000000
Delta_time = 1

# Computing VaR and ES with variance-covariance method
[VaR, ES] = AnalyticalNormalMeasures(alpha, w, notional, Delta_time, log_returns)
print("VaR of the ptf: ", VaR)
print("ES of the ptf:  ", ES)

# Plausibility check
VaR_plausibility = plausibilityCheck(log_returns, w, alpha, notional, 1)
print("VaR with plausibility check: ", VaR_plausibility)

# ----------------------------------------------------------------------------------------------------------------------------
# %%
# PUNTO 1

from utilities import *

# Importing data as above
data = pd.read_csv("EUROSTOXX50_Dataset.csv", index_col=0)

# Taking the stock prices for firms of the exercise 1-A
index = ["TTEF.PA", "AXAF.PA", "SASY.PA", "VOWG_p.DE"]
data1 = data[index].copy()
data1 = data1.fillna(method="ffill")
data1.index = pd.to_datetime(data1.index)

# Mining the data we need as above
current_date = datetime.datetime(2019, 3, 20)
d = 365 * 5 + 1
start_date = current_date - datetime.timedelta(d)
data1 = data1.loc[start_date:current_date]

# Computing log-return
log_returns = np.log(data1.shift(0) / data1.shift(1))
log_returns = log_returns.dropna()

# Getting last prices form the data
prices = data1.values
last_prices = prices[-1, :]

# Computing the weights
number_shares = np.array([25000, 20000, 20000, 10000])
w = number_shares * last_prices / np.dot(number_shares, last_prices)
alpha = 0.95
ptfValue = np.dot(last_prices, number_shares)

# Computing the VaR via Historical Simulation
[VaR, ES] = HSMeasurements(log_returns, alpha, w, np.dot(number_shares, last_prices), 1)
print("-----")
print("VaR of the ptf with HS: ", VaR)
print("ES of the ptf with HS:  ", ES)

# Plausibility check
print("VaR with plausibility: ", plausibilityCheck(log_returns, w, alpha, ptfValue, 1))

# Sampling random log returns
M = 200  # As M increases we obtain results closer and closer to the Historical Simulation method
bootstrap_returns = bootstrapStatistical(M, log_returns)

# Computing VaR with statistical bootstrap method (once we have sampled the log return we can simply call HSMeasurements)
[VaR_bootstrap, ES_bootstrap] = HSMeasurements(bootstrap_returns, 0.95, w, np.dot(number_shares, last_prices), 1)
print("-----")
print("VaR of the ptf with bootstrap: ", VaR_bootstrap)
print("ES of the ptf with bootstrap:  ", ES_bootstrap)

# Plausibility Check
print("VaR with plausibility: ", plausibilityCheck(bootstrap_returns, w, alpha, ptfValue, 1))

# %%
# WHS APPROACH
from utilities import *

# Importing data
data = pd.read_csv("EUROSTOXX50_Dataset.csv", index_col=0)

# Selecting firms of exercise 1_B
index = ["ADSGn.DE", "AIR.PA", "BBVA.MC", "BMWG.DE", "DTEGn.DE"]
data1 = data[index].copy()
data1 = data1.fillna(method="ffill")
data1.index = pd.to_datetime(data1.index)
data1 = data1.loc[start_date:current_date]

# Computing log-returns
log_returns = np.log(data1.shift(0) / data1.shift(1))
log_returns = log_returns.dropna()

# Creating the weights vector
w = np.array([0.20, 0.20, 0.20, 0.20, 0.20])

# Setting parameters for the WHS
alpha = 0.95
lamb = 0.95
ptfValue = 1

# Computing VaR and ES with WHS
[VaR_WHS, ES_WHS] = WHSMeasurements(log_returns, alpha, lamb, w, ptfValue, 1)
print("-----")
print("VaR of the ptf with WHS: ", VaR_WHS)
print("ES of the ptf with WHS:  ", ES_WHS)

# Plausability check
print("VaR with plausibility: ", plausibilityCheck(log_returns, w, alpha, ptfValue, 1))

#%%
# PUNTO 1.3

# Reading the file of indexes
name = pd.read_csv("_indexes.csv", index_col=0)

# Dropping the fourth & taking the first 18 companies
name = name.drop(4)
index = name.iloc[0:18, 0].values

# Importing data
data = pd.read_csv("EUROSTOXX50_Dataset.csv", index_col=0)
data1 = data[index].copy()
data1 = data1.fillna(method = "ffill")
data1.index = pd.to_datetime(data1.index)

# Computing the log-returns as above
log_returns = np.log(data1.shift(0)/data1.shift(1))
log_returns = log_returns.dropna()

# Computing the log returns yearly variance-covariance matrix and mean vector
var_cov_matrix = log_returns.cov().to_numpy()
yearlyCov = 256*var_cov_matrix
mu = log_returns.mean().to_numpy() #MEDIA
yearlyMu = 256*mu

# Computing the weights and setting parameters
weights = np.ones((1,18))*1/18
time_lag = 10/256
ptfValue = 1
alpha = 0.95

# Computing VaR and ES varying the number of components taken in consideration
print("--- Different values of VaR increasing the value of the number of principal components we consider ---")
print("")
for i in range(1, 6):
    VaR_PCA, ES_PCA = PrincCompAnalysis(yearlyCov, yearlyMu, weights, time_lag, alpha, i, ptfValue)
    print("VaR of the ptf with PCA: ", VaR_PCA)
    print("ES of the ptf with PCA:  ", ES_PCA)
    print("VaR with plausibility: ", plausibilityCheck(log_returns, weights, alpha, ptfValue, 10))
    print("-----")

# %%

# FULL MONTECARLO
# Importing data as above
data = pd.read_csv("EUROSTOXX50_Dataset.csv", index_col=0)
index = ["BMWG.DE"]
data1 = data[index].copy()
data1 = data1.fillna(method = "ffill")
data1.index = pd.to_datetime(data1.index)

# Selecting the time window we need from the stock prices
current_date = datetime.datetime(2017,1,16)
d = 365*2 + 1
start_date = current_date - datetime.timedelta(d);
expiry_date = datetime.datetime(2017,4,18)
data1 = data1.loc[start_date:current_date]  #data we need

# Computing the current stock price
prices = data1.values
last_price = prices[-1]

# Computing the number of stock in the portfolio
stockPosition = 1186680 # euro in BMW
number_stock = stockPosition/last_price #number of BMW in the ptf (equivalent to number of call we short)

# Setting parameters
strike = 25
vol = 0.154
dividend_yield = 0.031
risk_free_rate = 0.05
final_date = datetime.datetime(2017,4,18) # expiry_date
NumberOfDaysPerYears = 256
timeToMaturityInYears = (expiry_date-current_date).days/NumberOfDaysPerYears
riskMeasureTimeIntervalInYears = 10/NumberOfDaysPerYears
alpha = 0.95

# Computing log returns
log_returns = np.log(data1.shift(0)/data1.shift(1))
log_returns = log_returns.dropna()

print("-----")
# Computing the VaR with the full evaluation Monte Carlo we understand from the suggestions on the forum (number of iteration = 100k)
VaR_MC = FullMonteCarloVaR(log_returns, number_stock, number_stock, last_price, 25, risk_free_rate, dividend_yield, vol, timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha,NumberOfDaysPerYears)
print("Full evaluation Monte Carlo - VaR: ",VaR_MC)

#%%
# DELTA-NORMAL
data = pd.read_csv("EUROSTOXX50_Dataset.csv", index_col=0)
index = ["BMWG.DE"]
data1 = data[index].copy()
data1 = data1.fillna(method = "ffill")
data1.index = pd.to_datetime(data1.index)

# Selecting the time window we need from the stock prices
current_date = datetime.datetime(2017,1,16)
d = 365*2 + 1
start_date = current_date - datetime.timedelta(d);
expiry_date = datetime.datetime(2017,4,18)
data1 = data1.loc[start_date:current_date]  #data we need

# Computing the current stock price
prices = data1.values
last_price = prices[-1]

# Computing the number of stock in the portfolio
stockPosition = 1186680 # euro in BMW
number_stock = stockPosition/last_price #number of BMW in the ptf (equivalent to number of call we short)

# Setting parameters
strike = 25 #too low in our opinion to see some significant result
vol = 0.154
dividend_yield = 0.031
risk_free_rate = 0.05
expiry_date = datetime.datetime(2017,4,18) # expiry_date
NumberOfDaysPerYears = 256
timeToMaturityInYears = (expiry_date-current_date).days/NumberOfDaysPerYears
riskMeasureTimeIntervalInYears = 10/NumberOfDaysPerYears
alpha = 0.95

# Computing log returns
log_returns = np.log(data1.shift(0)/data1.shift(1))
log_returns = log_returns.dropna()

numberCall = number_stock

# Computing the portfolio VaR with Delta-Normal approach
VaR = DeltaNormalVaR(log_returns, number_stock, numberCall, last_price, strike, risk_free_rate, dividend_yield, vol, timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha, NumberOfDaysPerYears)
print("-----")
print("VaR with Delta-Normal method: ", VaR)
print("-----")

#%%

# CLIQUET OPTION
#Importing data regarding the bootstraped discount's curve (altready filtered in excel)
# "data_from_bootstrap" already contains discounts and zero rates up to seven years
bootstrap = pd.read_csv("data_from_bootstrap.csv", index_col=0)
bootstrap = bootstrap.dropna(axis=1) #dropping null on the column
# Index as datetime variables
bootstrap.index = pd.to_datetime(bootstrap.index)

survProbs_ISP = np.array([1.000000000000000, 0.995178319062266, 0.989368860449229, 0.982531958006141, 0.974039080328908, 0.966856037683655, 0.960934493204168, 0.952720422113471])

# Computing the defaultable discounts
default_B = bootstrap.values[:,0]*survProbs_ISP
#print(default_B)
# Computing forward discounts
nodef_B = bootstrap.values[:,0]
fwd_B = np.array([ nodef_B[i]/nodef_B[i-1] for i in range(1,len(nodef_B))])
#print(fwd_B)

# Computing yearfracs
# The vector yearFrac contains yearfrac between following years
yearFrac = yearfrac(bootstrap.index[0:-1],bootstrap.index[1:],6).to_numpy()
#print(yearFrac)
# The vector deltaT contains yearfrac between the settle and each year
deltaT = yearfrac(bootstrap.index[0],bootstrap.index[1:],3).to_numpy()

# Computing the forward rates
fwd_rate = - np.log(fwd_B) / yearFrac

# Setting parameters
sigma = 0.2
l = 0.99
k = 1/l
notional = 30000000
s0 = 1
recovery = 0.4

# Computing d1 and d2 according to the formula for the cliquet option
d1_bs = (math.log(l) + (fwd_rate + sigma**2/2)*yearFrac)/(sigma*np.sqrt(yearFrac)) #yearfrac???
d2_bs = d1_bs - sigma*np.sqrt(yearFrac)

# Computing the coupons' expected value
#print(fwd_rate)
#print(bootstrap.values[:,1])

### PERCHE FWD_RATE E NON I RATES DAL BOOTSTRAP
coupon = l*s0*np.exp(bootstrap.values[0,1]*deltaT)*(sc.stats.norm.cdf(d1_bs)-k*np.exp(-fwd_rate*yearFrac)*sc.stats.norm.cdf(d2_bs)) #cash flows #aggiunto un meno nell'exp

# Defining defaultProb that is usefull for the recovery term in the closed formula
defaultProb = np.array([survProbs_ISP[i-1] - survProbs_ISP[i] for i in range(1,len(survProbs_ISP))])

# Discounting the coupon
discountedCoupon = coupon*bootstrap.values[1:,0]

# The vector sum contains in the i-th the sum of the discounted coupons related to the following years
# It is used in the recovery term in the closed formula
discountedCoupon = np.flip(discountedCoupon)
sum = np.cumsum(discountedCoupon)
sum = np.flip(sum)

# Computing the Cliquet Price via closed formula
closePrice = np.dot(coupon,yearFrac*default_B[1:,]) + recovery*np.dot(defaultProb,sum)

price_byISP = np.sum(discountedCoupon)

print("The following results must be multiplied by the portfolio notional (€30 MLN): ")
print("Cliquet option price with closed formula: ",closePrice)
print("Cliquet option price with no counterparty risk (ISP): ",price_byISP)

numIter = 10000
MCPrice = cliquetOptionMonteCarlo(survProbs_ISP, bootstrap, s0, l, sigma, notional, recovery, numIter)
print("Cliquet option price with MC: ", MCPrice)