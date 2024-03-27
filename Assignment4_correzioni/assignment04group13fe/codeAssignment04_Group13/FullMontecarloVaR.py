import numpy as np
from blsprice import blsprice


def FullMontecarloVaR(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend, volatility,
                      timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha, NumberOfDaysPerYears, Nsim=1000000):

    """
    The function calculates the VaR for a nonlinear portfolio composed by a certain number of long shares and a certain
    number of calls, by using a full Monte Carlo VaR.
    the call option has the share as underlying asset.

    INPUTS:
        logReturns:                     the log returns of the shares. pandas dataframe with datetime dates as indexes
        numberOfShares:                 number of shares in the portfolio
        numberOfCalls:                  number of calls in the portfolio
        stockPrice:                     price of the stock at current date
        strike:                         strike price of the call option
        rate:                           risk-free rate
        dividend:                       dividend of the share
        volatility:                     volatility of the share
        timeToMaturityInYears:          time to maturity of the call option expressed in years
        riskMeasureTimeIntervalInYears: delta I'm interested in for the VaR
        alpha:                          significance level for VaR calculation
        numberOfDaysPerYears:           number of trading days in a year
        Nsim:                           number of simulations for the Monte Carlo approach

    OUTPUTS:
        VaR:                            the value at risk for the portfolio calculated for the risk measure time interval
    """

    lambda_ = 0.95 # Lambda used in the function that calculates the weights
    ticker = logReturns.columns[0]  # Get the first ticker name

    # Get the log returns in windows of 10 days
    logReturnsRollingWindow = logReturns.rolling(window=round(riskMeasureTimeIntervalInYears*NumberOfDaysPerYears)).sum().dropna()
    t = (logReturnsRollingWindow.index[-1] - logReturnsRollingWindow.index[0]).days
    historicalDates = np.array([(date - logReturnsRollingWindow.index[0]).days for date in logReturnsRollingWindow.index[:]])

    # Calculate the weights with a WHS approach
    C = (1 - lambda_) / (1 - lambda_ ** len(historicalDates))
    weights = np.array([C * lambda_ ** (t - s) for s in historicalDates])

    # Extract Nsim samples using the weights calculated above as pdf
    extractedS = np.random.choice(np.arange(0, len(historicalDates), 1), size=Nsim, p=weights / np.sum(weights))

    # Get the corresponding log returns as dX(t+delta)
    X = np.array(logReturnsRollingWindow[ticker].iloc[extractedS])

    # S(t+delta)
    prices = np.zeros((Nsim, 2))
    prices[:, 0] = stockPrice * np.exp(X)

    # C(t+delta)
    ttmInDelta = timeToMaturityInYears - riskMeasureTimeIntervalInYears  # used to get the price of the call at time t + delta
    prices[:, 1] = blsprice(prices[:, 0], strike, rate, dividend, volatility, ttmInDelta, flag=1)

    # value of the call that I'm shorting, today
    priceToday = blsprice(stockPrice, strike, rate, dividend, volatility, timeToMaturityInYears, flag=1)

    lossesSharesPortfolio = stockPrice - prices[:, 0]
    lossesDerivativePortfolio = priceToday - prices[:, 1]
    lossesTotalPortfolio = numberOfShares[0] * lossesSharesPortfolio + numberOfCalls * lossesDerivativePortfolio

    VaR = np.percentile(lossesTotalPortfolio, 100*(alpha))

    return VaR
