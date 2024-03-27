import numpy as np
from scipy.stats import norm

def blsprice(S, K, r, q, sigma, ttm, flag=1):
    """
    The function calculates the price of a European call or put option using the Black-Scholes model.

    INPUTS:
        S:      Current price of the underlying asset
        K:      Strike price of the option
        r:      Risk-free interest rate
        q:      Dividend yield of the underlying asset
        sigma:  Volatility of the underlying asset
        ttm:    Time to maturity of the option (in years)
        flag:   Type of the option (1 for call, -1 for put)

    OUTPUTS:
        price:  The price of the option
    """
    d_1 = (np.log(S/K) + (r-q+sigma**2/2)*ttm)/(sigma*np.sqrt(ttm))
    d_2 = d_1 - sigma*np.sqrt(ttm)
    return flag*S*np.exp(-q*ttm)*norm.cdf(flag*d_1) - flag*K*np.exp(-r*ttm)*norm.cdf(flag*d_2)

def blsdelta(stockPrice, strike, rate, dividend, volatility, timeToMaturityInYears):
    """
    The function calculates the delta of a European call option using the Black-Scholes model.

    INPUTS:
        stockPrice:             Current price of the underlying stock
        strike:                 Strike price of the option
        rate:                   Risk-free interest rate
        dividend:               Dividend yield of the underlying stock
        volatility:             Volatility of the underlying stock
        timeToMaturityInYears:  Time to maturity of the option (in years)

    OUTPUTS:
        delta:  The delta of the option
    """
    d1 = (np.log(stockPrice / strike) + (rate - dividend + 0.5 * volatility ** 2) * timeToMaturityInYears) / (volatility * np.sqrt(timeToMaturityInYears))
    delta = np.exp(-dividend * timeToMaturityInYears) * norm.cdf(d1)
    return delta
