# Functions - Assignment 4 Group 14

# Import the required packages
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.stats import t

def relevant_returns(stocks_data, indices_data, datestart, dateend, specific_names):
    # This function returns relevant stocks and their returns
    # Inputs:
    # stocks_data: Stocks data from the dataset
    # indices_data: Indexes from _indexes.csv
    # datestart: 20 o 23 Feb 2015
    # dateend: 20 Feb 2020
    # specific_names: Vector of the stocks considered

    # Compute relevant indices for time and companies

    # 'start' and 'end' are computed by finding the indices of datestart and dateend within the 'Date' column of stocks_data
    # Filtering the rows where the ‘Date’ column matches the string representation of datestart/dateend
    # .index[0] retrieves the index (row number) of the first matching row
    start = stocks_data.loc[stocks_data['Date'] == str(datestart)].index[0]
    end = stocks_data.loc[stocks_data['Date'] == str(dateend)].index[0]

    # Find the numeric index corresponding to each company name in specific_names
    # We add 1 to each index because Python uses 0-based indexing, while the stocks_indices should be 1-based (to match the actual column indices)
    stocks_indices = indices_data.loc[indices_data['Name'].isin(specific_names)].index + 1

    # Extract relevant stocks data using slicing
    # We selected rows from start to end (inclusive) and columns specified by stocks_indices
    # .reset_index(drop=True) resets the index of the resulting DataFrame
    relevant_stocks = stocks_data.iloc[start:end+1, stocks_indices].reset_index(drop=True)

    # Compute returns from the relevant stocks
    shifted_stocks = relevant_stocks.shift(1)
    returns = np.log(relevant_stocks[1:end + 1] / shifted_stocks[1:end +1])

    return [returns, relevant_stocks]

def AnalyticalNormalMeasures(alpha, weights, portfolioValue, riskMeasureTimeIntervalInDay, returns):
    # This function computes VaR and ES via t-Student parametric approach
    # Inputs:
    # alpha: Confidence level
    # weights: Weights of the assets in the portfolio
    # portfolioValue: Notional value of the portfolio
    # riskMeasureTimeIntervalInDay: 1
    # returns: returns of the stocks computed with function relevant_returns

    # Compute mean and covariance matrix of the returns
    mu = -returns.mean()  # mean of the loss
    cov = returns.cov()  # covariance
    delta = riskMeasureTimeIntervalInDay
    nu = 4  # Degrees of freedom of t-Student

    # Compute VaR
    VaR_std = t.ppf(alpha, nu)
    VaR = portfolioValue * (delta * np.dot(weights,mu) + np.sqrt(delta * np.dot(weights, np.dot(cov, weights))) * VaR_std)

    # Compute ES
    ES_std = (nu + (t.ppf(alpha, nu)) ** 2) / (nu-1) * 1 / (1 - alpha) * t.pdf(t.ppf(alpha, nu), nu)
    ES = portfolioValue * (delta * np.dot(weights,mu) + np.sqrt(delta * np.dot(weights, np.dot(cov, weights))) * ES_std)

    return [VaR, ES]

def Call_BLK_price(F, K, T, r, d, vol):
    d1 = (np.log(F/K) + (r-d+0.5*vol**2)*T) / np.sqrt(T*(vol**2))
    d2 = (np.log(F/K) + (r-d-0.5*vol**2)*T) / np.sqrt(T*(vol**2))
    price = F * np.exp(-d * T) * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
    return price

# Bootstrap method
# Creating a function for samples
def bootstrapStatistical(numberOfSamplesToBootstrap, returns):
    # Generate random indices
    indices = np.random.choice(np.arange(returns.shape[0]), size=numberOfSamplesToBootstrap, replace=True)
    # Create an array of samples using the selected indices
    samples = returns[indices, :]

    return samples

def HSMeasurements(returns, alpha, weights, portfolioValue, riskMeasureTimeIntervalInDay):
    # This function returns the Historical simulation calculation Var and ES
    # Inputs:
    # returns: Stocks returns
    # alpha: quantile
    # weights: weights of the stocks
    # portfolioValue: Value of the portfolio now
    # riskMeasureTimeIntervalInDay: time interval

    # Computing the losses under the empirical set of realized values for factors
    L = []
    n = len(returns)
    for i in range(len(returns)):
        L.append(-portfolioValue * np.dot(weights, returns[i, :]))

    L.sort(reverse=True)
    VaR = L[int(n * (1 - alpha))-1 ]
    ES = np.mean(L[0:int(n * (1 - alpha))])

    return VaR, ES

def WHSMeasurements(returns, alpha, lam, weights, portfolioValue,riskMeasureTimeIntervalInDay):
    # This function returns the Weighted Historical simulation calculation Var and ES
    # Inputs:
    # returns: Stocks returns
    # alpha: quantile
    # lambda: in order to give less importance to past losses
    # weights: weights of the stocks
    # portfolioValue: Value of the portfolio now
    # riskMeasureTimeIntervalInDay: time interval

    C = (1 - lam) / (1 - lam ** len(returns))

    # Calculate weights
    w = [C * lam ** (len(returns) - 1 - i) for i in range(len(returns))]
    w = np.array(w)

    # Initialize L
    L =[]

    for i in range(len(returns)):
        L.append(-portfolioValue * np.dot(weights, returns[i, :]))

    L = np.array(L)

    # Sort in descending order based on L
    w = w[np.argsort(L)[::-1]]
    L = L[np.argsort(L)[::-1]]

    # Initialize variables
    cumulative_sum = 0
    index = 0

    # Accumulate values until cumulative sum reaches 1 - alpha
    while cumulative_sum <= 1 - alpha:
        cumulative_sum += w[index]
        index += 1

    VaR = L[index-1]
    ES = np.dot(w[:index-1],L[:index-1])/(cumulative_sum-w[index-1])

    return VaR, ES


def PrincCompAnalysis(yearlyCovariance, yearlyMeanReturns, weights, H, alpha,numberOfPrincipalComponents, portfolioValue):
    # PCA for VaR and ES
    # yearlyCovariance: covariance matrix
    # yearlyMeanReturns
    # weights: weights of each stock in the portfolio
    # H: delta in years
    # alpha: quantile VaR
    # numberOfPrincipalComponents
    # portfolioValue

    # Sorting the values in descending ordre respect to the eigenvalues
    eigenval, eigenvect=np.linalg.eig(yearlyCovariance)
    eigenvect = eigenvect[:, np.argsort(eigenval)[::-1]].T
    eigenval=eigenval[np.argsort(eigenval)][::-1]


    # Computing mu_hat and sigma_hat
    mu_hat=np.dot(eigenvect, yearlyMeanReturns)
    weights_hat=np.dot(eigenvect, weights)

    # Computing the reduced mu and sigma^2 of the portfolio
    mu_red=np.dot(weights_hat[:numberOfPrincipalComponents],mu_hat[:numberOfPrincipalComponents])
    sigma2_red=np.dot(weights_hat[:numberOfPrincipalComponents]**2,eigenval[:numberOfPrincipalComponents])

    # We compute the VaR(mu+sqrt(H)*sigma*VaR_std) and ES(mu+sqrt(H)*sigma*ES_std)
    VaR=portfolioValue*(-mu_red*H+np.sqrt(H*sigma2_red)*norm.ppf(alpha))
    ES=portfolioValue*(-mu_red*H+np.sqrt(H*sigma2_red)*(1/(1-alpha))*(norm.pdf(norm.ppf(alpha))))
    return VaR,ES

def plausibilityCheck(returns, portfolioWeights, alpha, portfolioValue, riskMeasureTimeIntervalInDay):
    # sensitivity=weights because of linearity

    C = np.corrcoef(returns.T)
    l = np.quantile(-returns, (1 - alpha) , axis=0)
    u = np.quantile(-returns, alpha , axis=0)
    sVaR=portfolioWeights*((abs(l)+abs(u))/2)

    VaR = portfolioValue * np.sqrt(np.dot(np.dot(sVaR, C), sVaR)) * np.sqrt(riskMeasureTimeIntervalInDay)
    return VaR

def FullMonteCarloVaR(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend, volatility,
                      timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha, NumberOfDaysPerYears, l, niter):
    # Inputs:
    # logReturns: Logarithmic Returns over VaR time Horizon
    # numberOfShares: Number of shares held in portfolio
    # numberOfCalls: Number of calls held in portfolio
    # stockPrice: Stock initial price
    # strike: Option strike price
    # rate: Risk Free rate
    # dividend: Option continuous dividend yield
    # volatility: Option volatility
    # timeToMaturityInYears: Option TTM
    # riskMeasureTimeIntervalInYears: VaR Time horizon in years
    # alpha:
    # NumberOfDaysPerYears: Self explanatory
    # l : Lambda Historical weighted simulation parameter
    # niter: Number of historical data sampled
    # seed: Seed for consistent sampling

    # Vector containing integers from 1 to total number of log returns
    delta_time = np.arange(1, len(logReturns) + 1)

    # Define weights vector
    weights = np.array([(1 - l) / (1 - l ** len(logReturns)) * l ** (len(logReturns) - delta) for delta in delta_time])

    # Sample niter returns and weights from uniform distribution in a vector
    indexes = np.random.randint(0, len(logReturns), int(niter))
    Sampled_Returns = np.array(logReturns.iloc[indexes])

    # Matching weights to log returns
    Weights = weights[indexes]

    # Simulate niter future stock and call prices at the end of VaR time horizon based on niter historical data sampled
    StockPriceNdays = (stockPrice * np.exp(Sampled_Returns))

    Call_today = Call_BLK_price(stockPrice, strike, timeToMaturityInYears, rate, dividend, volatility)
    CallPriceNdays = []

    for i in range(niter):
        CallPriceNdays = np.append(CallPriceNdays, Call_BLK_price(StockPriceNdays[i], strike,
                                                    timeToMaturityInYears - riskMeasureTimeIntervalInYears,
                                                      rate, dividend, volatility))

    # Compute the total loss in niter different scenarios as sum of Stock losses and Call losses
    Loss_S = -numberOfShares * (StockPriceNdays.T - stockPrice)
    LossCall = -numberOfCalls * (CallPriceNdays - Call_today)
    TotalLoss = Loss_S + LossCall

    # Sort the loss vector and corresponding weight vector
    L_sorted = -np.sort(-TotalLoss)

    indexes2 = np.argsort(TotalLoss)[::-1]
    weights_sorted = Weights[indexes2]

    # Compute 1-alpha quantile
    index = np.cumsum(weights_sorted) <= (1 - alpha) * np.sum(weights_sorted)

    index_star = np.sum(index)

    VaR = L_sorted[:, index_star - 1]

    return VaR


def DeltaNormalVaR(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend, volatility,
                   timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha, NumberOfDaysPerYears, l, niter):
    # Inputs:
    # logReturns: Logarithmic Returns over VaR time Horizon
    # numberOfShares: Number of shares held in portfolio
    # numberOfCalls: Number of calls held in portfolio
    # stockPrice: Stock initial price
    # strike: Option strike price
    # rate: Risk Free rate
    # dividend: Option continuous dividend yield
    # volatility: Option volatility
    # timeToMaturityInYears: Option TTM
    # riskMeasureTimeIntervalInYears: not used
    # alpha:
    # NumberOfDaysPerYears
    # l: Lambda Historical weighted simulation parameter
    # niter: Number of historical data sampled

    # Compute weights
    delta = np.arange(1, len(logReturns) + 1)

    # Define weights vector
    weights = np.array([(1 - l) / (1 - l ** len(logReturns)) * l ** (len(logReturns) - delta) for delta in delta])

    # Sample niter returns and weights from uniform distribution in a vector
    indexes = np.random.randint(0, len(logReturns), int(niter))

    Sampled_Returns = np.array(logReturns.iloc[indexes])

    # Match weights to log returns
    Weights = weights[indexes]

    # Simulate Stock Price at the end of VaR time horizon and the Loss Registered
    StockPriceNdays = stockPrice * np.exp(Sampled_Returns)
    Loss_S = -numberOfShares * (StockPriceNdays - stockPrice)

    # Approximate Call Price change by first order taylor expansion
    delta = Call_BLK_price(stockPrice, strike, timeToMaturityInYears, rate, dividend, volatility)
    Loss_Call = -numberOfCalls * delta * (StockPriceNdays - stockPrice)

    # Compute total loss
    TotalLoss = Loss_S + Loss_Call

    # Sort Loss Vector and matching weights vector
    L_sorted = -np.sort(-TotalLoss.flatten())
    indexes2 = np.argsort(TotalLoss.flatten())[::-1]
    weights_sorted = Weights[indexes2]

    # Compute quantile
    index = np.cumsum(weights_sorted) <= (1 - alpha) * np.sum(weights_sorted)
    index_star = np.sum(index)
    VaR = L_sorted[index_star - 1]

    return VaR

def cliquet_price_MC(settlement_date, payoff_dates, notional, L, S0, volatility, recovery, N):
    # Compute the price of a cliquet option given: notional, yearly payoff dates (payoff =  [𝐿 ∗ 𝑆(𝑡_i) − 𝑆(𝑡_i-1)]+),
    # where L is the participation coefficient and S is the underlying equity with given constant volatility
    # Using MC methods with antithetic variables

    # INPUT
    # settlement
    # payoff_dates
    # notional
    # L:                    Participation coefficient.
    # SO:                   Underlying initial value
    # volatility
    # recovery
    # N:                    Number of simulation of the underlying dynamic

    # OUTPUT
    # price_rf                Simulated price of the option (risk-free)
    # price_def               Simulated price of the option (considering the risk)

    # Read data
    bootstrap_discounts = pd.read_csv('Discounts.csv').to_numpy()
    bootstrap_survprobs = pd.read_csv('ISP_survival_probabilities.csv').to_numpy()

    # Concatenate settlement date with payoff dates
    dates = np.concatenate(([settlement_date], payoff_dates))

    # Calculate the difference in days between consecutive dates and convert them to fraction using 30/360
    yearfrac = (np.diff(dates).astype('timedelta64[D]').astype(int)) / 360

    # Compute relevant discounts, rates, survival probabilities and forward discounts
    discounts = np.interp(payoff_dates, bootstrap_discounts[:, 0], bootstrap_discounts[:, 1])
    discounts = np.insert(discounts, 0, 1)
    fw_rates = -np.log(discounts[1:] / discounts[:-1]) / yearfrac
    survprobs = np.interp(payoff_dates, bootstrap_survprobs[:, 0], bootstrap_survprobs[:, 1])
    def_discounts = discounts[1:]*survprobs

    # Underlying S dynamic
    # For year 0 all the price are equal to S0
    S = np.zeros((len(dates), 2*N))
    S[0, :] = S0*np.ones(2*N)
    # For year from 1 to 7 we simulate iteratively the values of the underlying using the GBM dynamics
    for i in range(1, len(dates)):
        # Use the antithetic variables technique to speed up the process
        g = np.random.randn(N)
        g = np.concatenate((g, -g))
        S[i, :] = S[i-1, :]*np.exp((fw_rates[i-1]-(volatility**2)/2) * yearfrac[i-1] + np.sqrt(yearfrac[i-1]) * volatility * g)
    # Compute the yearly payoffs with the given formula for each simulation
    payoff = np.maximum(L*S[1:, :] - S[0:-1, :], 0)
    # Compute the yearly mean of all the payoffs simulated
    payoff = np.mean(payoff, axis=1)

    # Compute the "risk-free" price
    price_rf = np.sum(discounts[1:]*payoff)

    # Compute the "risky" price
    price_def = np.sum(def_discounts*payoff*(survprobs + recovery*(1 - survprobs)))

    return [price_rf, price_def]


def cliquet_price_analytic(settlement_date, payoff_dates, notional, L, volatility, recovery):
    # Compute the price of a cliquet option given: notional, yearly payoff dates (payoff =  [𝐿 ∗ 𝑆(𝑡_i) − 𝑆(𝑡_i-1)]+),
    # where L is the participation coefficient and S is the underlying equity with given constant volatility
    # Using analytical formula

    # INPUT
    # settlement
    # payoff_dates
    # notional
    # L:                    Participation coefficient.
    # volatility
    # recovery

    # OUTPUT
    # price_rf              Exact price of the option (risk-free)
    # price_def             Exact price of the option (considering the risk)

    # Read data
    bootstrap_discounts = pd.read_csv('Discounts.csv').to_numpy()
    bootstrap_survprobs = pd.read_csv('ISP_survival_probabilities.csv').to_numpy()

    # Concatenate settlement date with payoff dates
    dates = np.concatenate(([settlement_date], payoff_dates))

    # Calculate the difference in days between consecutive dates and convert them to fraction using 30/360
    yearfrac = (np.diff(dates).astype('timedelta64[D]').astype(int))/360

    # Compute relevant discounts, rates, survival probabilities and forward discounts
    discounts = np.interp(payoff_dates, bootstrap_discounts[:, 0], bootstrap_discounts[:, 1])
    discounts = np.insert(discounts, 0, 1)
    fw_rates = -np.log(discounts[1:]/discounts[:-1])/yearfrac
    rates = -np.log(discounts[1:])/yearfrac
    fw_discounts = discounts[1:]/discounts[:-1]
    survprobs = np.interp(payoff_dates, bootstrap_survprobs[:, 0], bootstrap_survprobs[:, 1])

    # Compute the EU Call prices using Black76 formula
    call_prices = np.zeros(len(payoff_dates))
    for i in range(0, len(call_prices)):
        call_prices[i] = Call_BLK_price(fw_discounts[i], 1/L, yearfrac[i], fw_rates[i],0, volatility)

    # Compute the "risk-free" price
    price_rf = L*sum((1/discounts[0:-1])*call_prices)

    # Compute the "risky" price
    price_def = L*sum((1/discounts[0:-1])*call_prices*(survprobs + recovery*(1 - survprobs)))

    return [price_rf, price_def]