
import numpy as np
import pandas as pd

# function to compute the historical simulation VaR and ES

def HSMeasurements(returns, alpha, weights, portfolioValue, riskMeasureTimeIntervalInDay):

    """
        This function computes the historical simulation VaR and ES for a given portfolio of assets.

        Args:
        - returns: a pandas dataframe of asset returns (n, m)
        - alpha: the confidence level (scalar)
        - weights: a numpy array of portfolio weights (m, )
        - portfolioValue: the value of the whole portfolio (scalar)
        - riskMeasureTimeIntervalInDay: the time interval for risk measure expressed in days (scalar)

    """
    
    losses = - portfolioValue * (returns @ weights)
    losses = losses.sort_values(ascending=False)

    n_quantile = int( (1-alpha) * len(losses) )

    VaR = losses.iloc[n_quantile] * np.sqrt(riskMeasureTimeIntervalInDay)

    ES = losses.iloc[:n_quantile].mean() * np.sqrt(riskMeasureTimeIntervalInDay)

    return ES, VaR

def bootstrapStatistical(numberOfSamplesToBootstrap, returns):
    """
        This function computes the bootstrap samples of the returns.
        
        Args:
        - numberOfSamplesToBootstrap: the number of samples to draw from the returns (integer)
        - returns: a pandas dataframe of asset returns (n, m)
    """

    return returns.sample(n=numberOfSamplesToBootstrap, replace=True)

def plausibilityCheck(returns, portfolioWeights, alpha, portfolioValue, riskMeasureTimeIntervalInDay):

    """
        This function computes the plausibility check for the VaR of a portfolio of assets.

        Args:
        - returns: a pandas dataframe of asset returns (n, m)
        - portfolioWeights: a numpy array of portfolio weights (m, )
        - alpha: the confidence level (scalar)
        - portfolioValue: the value of the whole portfolio (scalar)
        - riskMeasureTimeIntervalInDay: the time interval for risk measure expressed in days (scalar)
    """

    # extract the two quantiles
    lb = returns.quantile(1-alpha)
    ub = returns.quantile(alpha)

    # compute the signed VaR
    sVaR = portfolioValue * portfolioWeights * (abs(lb) + abs(ub)) / 2

    # correlation matrix
    corr = returns.corr()

    # compute the VaR with the thumb rule
    VaR = np.sqrt(riskMeasureTimeIntervalInDay) * np.sqrt(sVaR @ corr @ sVaR)

    return VaR