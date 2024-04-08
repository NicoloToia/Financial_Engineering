import numpy as np
import pandas as pd
import scipy.io

from Portfolio import Portfolio
import tickersGetter as tg
from CliquetOption import CliquetOption
from datetime import datetime
import matplotlib.pyplot as plt

# !!!: global variables are usually not a good idea, should really be a class attribute
global numberOfDaysPerYears
numberOfDaysPerYears = 256

np.random.seed(0)

def newExercise(name):
    print("\n\n\n\n=============================================================")
    print(f"Exercise {name}")

tickers = pd.read_csv('EUROSTOXX50/_indexes.csv', index_col=0)
prices = pd.read_csv('EUROSTOXX50/EUROSTOXX50_Dataset.csv', index_col=0)

# Fill missing data
prices.ffill(inplace=True)
# TODO: ERROR! No backward fill
prices.bfill(inplace=True)
# Convert the string index in dates
prices.index = pd.to_datetime(prices.index, format='%Y-%m-%d')


# Exercise 0
newExercise("0")
assetsZero = ['Adidas', 'Allianz', 'Munich Re', "L'Oréal"]
assetVec = tg.getTickers(tickers, assetsZero)
portfolio0 = Portfolio(assetVec, prices, np.array([0.25, 0.25, 0.25, 0.25]), 15000000)
# !!!: second date is a non-trading day
portfolio0.setDates(pd.to_datetime('2015-02-20'), pd.to_datetime('2020-02-20'))

VaR0 = portfolio0.VaR(0.99, 1, "t", 4)
ES0 = portfolio0.ES(0.99, 1, "t", 4)
portfolio0.printPortfolioDetails('t Student', VaR0, ES0)

# Exercise 1
# 1.A
newExercise("1.A")
assetsA = ['Total', 'AXA', 'Sanofi', 'Volkswagen']
assetVecA = tg.getTickers(tickers, assetsA)
n_shares = [25000, 20000, 20000, 10000]
portfolioA = Portfolio(assetVec, prices, n_shares)
portfolioA.setDates(pd.to_datetime('2014-03-20'), pd.to_datetime('2019-03-20'))

VaRA, ESA = portfolioA.HSM(0.95, 1, True)
VarABS = portfolioA.bootstrapStat(200, 0.95, 1, print=True)

VaRApcheck = portfolioA.pCheck(0.95, 1, print=True)

# 1.B
newExercise("1.B")
assetsB = ['Adidas', 'Airbus', 'BBVA', 'BMW', 'Deutsche Telekom']
assetVecB = tg.getTickers(tickers, assetsB)
weightsB = np.full(5, 1/5)
portfolioB = Portfolio(assetVecB, prices, weightsB, 1)
portfolioB.setDates(pd.to_datetime('2014-03-20'), pd.to_datetime('2019-03-20'))
VarB, ESB = portfolioB.WHSM(0.95, 0.95, 1, print=True)
VarBpcheck = portfolioB.pCheck(0.95, 1, print=True)

# 1.C
newExercise("1.C")
tick = tickers['Ticker']
assetVecC = (np.concatenate((tick[0:3], tick[4:19]))).tolist()
weightsC = np.full(18, 1/18)
portfolioC = Portfolio(assetVecC, prices, weightsC, 1)
portfolioC.setDates(pd.to_datetime('2014-03-20'), pd.to_datetime('2019-03-20'))

VaRC = portfolioC.pCheck(0.95, 10, print=True)

print("\nn | VaR (95%)        | ES (95%)")
print("--|------------------|-----------------")
for n in range(1, 6):
    var, es = portfolioC.PCA(0.95, n, 10/256)  # Calcolo del VaR e ES usando PCA
    print(f"{n:<2}| € {var:<15.6f}| € {es:<15.6f}")


# 2.A-B
newExercise("2.A-B")
assets2 = ['BMW']
endDate2 = pd.to_datetime('2017-01-16')
expiryDate = pd.to_datetime('2017-04-18')
assetVec2 = tg.getTickers(tickers, assets2)
portfolio2Value = 1186680
portfolio2 = Portfolio(assetVec2, prices, np.array([1]), portfolio2Value)
portfolio2.setDates(prices.index[0], endDate2)
portfolio2.updateNShares()
riskMeasureTimeIntervalInYears = 10/numberOfDaysPerYears
portfolio2.addEuOption(-portfolio2.nShares[0], 'call', assetVec2,  25, 0.005, 0.031, expiryDate, 0.154)

VaR2a = portfolio2.FullMC('EU_call_BMWG.DE', riskMeasureTimeIntervalInYears, 0.95, numberOfDaysPerYears, print=True)
VaR2b = portfolio2.DNVaR('EU_call_BMWG.DE',riskMeasureTimeIntervalInYears, 0.95, numberOfDaysPerYears, print=True)

# 3
newExercise("3")
print("-------------------------------------------------------------")
survProb = scipy.io.loadmat('survProb.mat')  # survival probabilities of ISP for years 1, ..., 7
DiscountCurve = scipy.io.loadmat('DiscountCurve.mat')  # discount curve with dates and discount factors
paymentDates = scipy.io.loadmat('paymentDates.mat')  # payment dates of payoffs
yearlySurvProb = survProb['survProb']
DiscountCurve = DiscountCurve['DiscountCurve']
dates = paymentDates['paymentDates']
L, volatility, notional, maturity = 0.99, 0.2, 30e6, 7,
today = DiscountCurve[0, 0]
recovery = 0.4

cliquetOptionPrice = CliquetOption(yearlySurvProb, dates, L, notional, today, maturity, DiscountCurve, volatility, recovery)
print(f"7y Cliquet Option price: {cliquetOptionPrice:>15.4f} €")

plt.show()