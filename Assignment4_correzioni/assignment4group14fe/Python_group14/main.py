# MAIN ASSIGNMENT 4 - GROUP 14

# Import the required packages
import numpy as np
import pandas as pd
import datetime
import functions

########################################################################################################################

# Set the seed for reproducibility
np.random.seed(14)

# Import data
stocks_data = pd.read_csv('EUROSTOXX50_Dataset.csv')
stocks_data.fillna(axis=0, method='ffill', inplace=True) # setting missing values to the previous available one
indices_data = pd.read_csv('_indexes.csv')

########################################################################################################################
# Exercise 0

# Parameters
datestart = datetime.date(2015, 2, 20) # prendo 20 o 23 Feb 2015?
dateend = datetime.date(2020, 2, 20)
specific_names = ['Adidas', 'Allianz', 'Munich Re', 'L\'Oréal']
alpha = 0.99 # Confidence Level
portfolioValue = 15e6 # Notional value
riskMeasureTimeIntervalInDay = 1

# Compute weights and returns
weights = np.ones(len(specific_names))/len(specific_names)
returns = functions.relevant_returns(stocks_data, indices_data, datestart, dateend, specific_names)[0]

# Compute VaR and ES
[VaR_1, ES] = functions.AnalyticalNormalMeasures(alpha, weights, portfolioValue, riskMeasureTimeIntervalInDay, returns)
print('Exercise 0 \n VaR: ', VaR_1, ', ES:', ES, '\n')

########################################################################################################################
# Exercise 1
#HS Simulation

# Initialisation
datestart = datetime.date(2014, 3, 20)
dateend = datetime.date(2019, 3, 20)
specific_names = ['TotalEnergies', 'AXA', 'Sanofi', 'Volkswagen Group']
alpha=0.95

# Calculate the returns
returns = functions.relevant_returns(stocks_data, indices_data, datestart, dateend, specific_names)[0]
shares = functions.relevant_returns(stocks_data, indices_data, datestart, dateend, specific_names)[1]
shares=shares.to_numpy()
returns=returns.to_numpy()
shares_t=shares[-1,:] #values shares now
n_shares=np.array([25000,20000, 20000,10000])
portfolioVal=np.dot(shares_t,n_shares)
w=(shares_t*n_shares)/portfolioVal
[VaR,ES]=functions.HSMeasurements(returns, alpha, w, portfolioVal, riskMeasureTimeIntervalInDay)
print('Exercise 1.1 \n ')
print('HS VaR: ', VaR, ', ES:', ES, '\n')
# We do the Bootsrap to get 200 samples
numberOfSamplesToBootstrap=200
returns_samples=functions.bootstrapStatistical(numberOfSamplesToBootstrap, returns)
[VaR,ES]=functions.HSMeasurements(returns_samples, alpha, w, portfolioVal, riskMeasureTimeIntervalInDay)
print('HS Bootstrap VaR: ', VaR, ', ES:', ES, '\n')
VaR_check = functions.plausibilityCheck(returns, w, alpha, portfolioVal, riskMeasureTimeIntervalInDay)
print('VaR Plausability check: ', VaR_check, '\n')

# WHS Simulation
# Initialisation
lam=0.95
alpha=0.95
WHSportfolioValue=1

specific_names = ['Adidas', 'Airbus', 'BBVA', 'BMW','Deutsche Telekom']
w = 1/5*np.ones(5)  #equally weighted equity
returns = functions.relevant_returns(stocks_data, indices_data, datestart, dateend, specific_names)[0]
returns = returns.to_numpy()
print('Exercise 1.2 \n ')

[VaR,ES] = functions.WHSMeasurements(returns, alpha, lam, w, WHSportfolioValue,riskMeasureTimeIntervalInDay)
print('WHS VaR: ', VaR, ', ES:', ES, '\n')

VaR_check = functions.plausibilityCheck(returns, w, alpha, WHSportfolioValue, riskMeasureTimeIntervalInDay)
print('VaR Plausability check: ', VaR_check, '\n')

#PCA
#Initialisation
PCAportfolioValue = 1  # value of the portfolio
H = 10  # time interval in days
n_companies = 18
alpha=0.95

names_PCA = indices_data['Name'][:n_companies+1].tolist()
names_PCA.remove('Adyen')  #we remove it because it has no data
w = 1/18*np.ones(18)  #equally weighted equity
returns = functions.relevant_returns(stocks_data, indices_data, datestart, dateend, names_PCA)[0]
returns = returns.to_numpy()
covariance = np.cov(returns.T)
meanreturns = np.mean(returns.T,1)
numberOfPrincipalComponents = 5
print('Exercise 1.3 PCA \n ')

for n in range(1,numberOfPrincipalComponents+1):
    [VaR, ES] = functions.PrincCompAnalysis(covariance, meanreturns, w, H, alpha,n, PCAportfolioValue)
    print('number of components:', n, 'VaR: ', VaR, ', ES:', ES, '\n')

VaR_check = functions.plausibilityCheck(returns, w, alpha, PCAportfolioValue, H)
print('Excercise 1.3 PCA, VaR Plausability check: ', VaR_check, '\n')
########################################################################################################################

# Exercise 2

Notional = 1186680
specific_name = ['BMW']
datecall = datetime.date(2017, 4, 18)
dateend = datetime.date(2017, 1, 16)
datestart = datetime.date(2015, 1, 16)
datetomorrow = datetime.date(2017, 1, 17)
strike = 25
dividend = 3.1/100
volatility = 15.4/100
NumberOfDaysPerYears = 365
riskMeasureTimeIntervalInYears = 10/NumberOfDaysPerYears
start = stocks_data.loc[stocks_data['Date'] == str(datestart)].index[0]
end = stocks_data.loc[stocks_data['Date'] == str(dateend)].index[0]
timeToMaturityInYears = (datecall - datetomorrow).days/NumberOfDaysPerYears
alpha = 0.95
rate = 0.5/100
lam = 0.95
niter = 10000

# Compute returns
[returns, stocks] = functions.relevant_returns(stocks_data, indices_data, datestart, dateend, specific_name)
stockPrice = stocks.iloc[13]["BMWG.DE"]  #position of BMW index

numberOfShares = int(Notional/stockPrice)
numberOfCalls = numberOfShares

VaR_2 = functions.FullMonteCarloVaR(returns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend, volatility, timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha,NumberOfDaysPerYears,lam,niter)
print('Exercise 2 \n Var with Full MonteCarlo: ', VaR_2)

VaR_3 = functions.DeltaNormalVaR(returns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend, volatility, timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha, NumberOfDaysPerYears, lam, niter)
print(' Var with DeltaNormal: ', VaR_3)

########################################################################################################################

# Exercise 3

# Parameters
settlement_date = 733457
payoff_dates = [733823, 734188, 734555, 734919, 735284, 735649, 736014]
notional = 30*(10**6)
L = 0.99
recovery = 0.4
S0 = 1
N = 10**6
volatility = 0.2

print('\nExercise 3: \n ')

# Compute the Cliquet price via closed formula
[price_rf_exact, price_def_exact] = functions.cliquet_price_analytic(settlement_date, payoff_dates, notional, L, volatility, recovery)
print("Risk-free exact price = ", price_rf_exact)
print("Risky exact price = ", price_def_exact)

[price_rf_MC, price_def_MC] = functions.cliquet_price_MC(settlement_date, payoff_dates, notional, L, S0, volatility, recovery, N)
print("Risk-free MC price = ", price_rf_MC)
print("Risky MC price = ", price_def_MC)