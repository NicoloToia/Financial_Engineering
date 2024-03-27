import numpy as np
from WHSMeasurements import WHSMeasurements
import blsprice as bls

def DeltaNormalVaR(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend,
volatility, timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha, NumberOfDaysPerYears, flag=1, Nsim=1000000):
    """
    The function calculates the Value at Risk (VaR) for a portfolio composed of a certain number of shares and call options,
    using either an analytical Weighted Historical Simulation (WHS) approach or a Monte Carlo approach.

    INPUTS:
        logReturns:                          the log returns of the shares. pandas dataframe with datetime dates as indexes
        numberOfShares:                      number of shares in the portfolio
        numberOfCalls:                       number of calls in the portfolio
        stockPrice:                          price of the stock at current date
        strike:                              strike price of the call option
        rate:                                risk-free rate
        dividend:                            dividend of the share
        volatility:                          volatility of the share
        timeToMaturityInYears:               time to maturity of the call option expressed in years
        riskMeasureTimeIntervalInYears:      delta I'm interested in for the VaR
        alpha:                               significance level for VaR calculation
        NumberOfDaysPerYears:                number of trading days in a year
        flag:                                flag to choose between analytical WHS approach (0) and Monte Carlo approach (1)
        Nsim:                                number of simulations for the Monte Carlo approach

    OUTPUTS:
        VaR:                                 the value at risk for the portfolio calculated for the risk measure time interval
    """

    # calculate the delta of the call in t_0
    deltaCall = bls.blsdelta(stockPrice, strike, rate, dividend, volatility, timeToMaturityInYears)
    portfolioValueSharesOnly = numberOfShares * stockPrice
    lambda_ = 0.95  # decay factor for WHS
    if flag == 0:  # analytical WHS approach
        VaR, ES = WHSMeasurements(logReturns, alpha, lambda_, np.array([1]), portfolioValueSharesOnly, riskMeasureTimeIntervalInYears*NumberOfDaysPerYears)

        DeltaNormalVaR = VaR * (1 + deltaCall*numberOfCalls/numberOfShares)

        return DeltaNormalVaR[0]

    else:  # MonteCarlo approach
        lambda_ = 0.95  # Lambda used in the function that calculates the weights
        ticker = logReturns.columns[0]  # Get the first ticker name

        # Get the log returns in windows of 10 days
        logReturnsRollingWindow = logReturns.rolling(window=10).sum().dropna()
        t = (logReturnsRollingWindow.index[-1] - logReturnsRollingWindow.index[0]).days
        historicalDates = np.array(
            [(date - logReturnsRollingWindow.index[0]).days for date in logReturnsRollingWindow.index[:]])

        # Calculate the weights with a WHS approach
        C = (1 - lambda_) / (1 - lambda_ ** len(historicalDates))
        weights = np.array([C * lambda_ ** (t - s) for s in historicalDates])

        # Extract Nsim samples using the weights calculated above as pdf
        extractedS = np.random.choice(np.arange(0, len(historicalDates), 1), size=Nsim, p=weights / np.sum(weights))

        # Get the corresponding log returns as dX(t+delta)
        X = np.array(logReturnsRollingWindow[ticker].iloc[extractedS])

        # S(t+delta)
        prices = stockPrice * np.exp(X)

        lossesSharesPortfolio = stockPrice - prices
        lossesTotalPortfolio = lossesSharesPortfolio*(1*numberOfShares + deltaCall*numberOfCalls)
        VaR = np.percentile(lossesTotalPortfolio, 100 * (alpha))
        return VaR
