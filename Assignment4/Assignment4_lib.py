
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
    returns_new = []
    for i in range(1, numberOfSamplesToBootstrap):
        # extract a random number in range from 1 to the number of rows of the returns
        n = np.random.randint(1, len(returns))
        returns_new = returns_new + returns.sample(n, replace=True)
    
    return returns_new

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

def WHSMeasurements(returns, alpha, lmb, weights, portfolioValue, riskMeasureTimeIntervalInDay):

    """
        This function computes the weighted historical simulation VaR and ES for a given portfolio of assets.

        Args:
        - returns: a pandas dataframe of asset returns (n, m)
        - alpha: the confidence level (scalar)
        - lmb = lambda
        - weights: a numpy array of portfolio weights (m, )
        - portfolioValue: the value of the whole portfolio (scalar)
        - riskMeasureTimeIntervalInDay: the time interval for risk measure expressed in days (scalar)

    """
    
    losses = - portfolioValue * (returns @ weights)
    losses = losses.sort_values(ascending=False)
    C = (1-lmb)/(1-lmb^(len(losses)))
    n = len(losses)
    sumWeights = 0
    weight = []
    for jj in range(n-1): 
        w = C* lmb^(n-jj)
        sumWeights = sumWeights+w
        if sumWeights>1-alpha:
            k = jj-1
            break
        else :
            weight[jj] = w
                
    VaR = losses.iloc[k] * np.sqrt(riskMeasureTimeIntervalInDay)

    ES = sum(losses.iloc[:k] @ weights)/sum(weights[:k]) * np.sqrt(riskMeasureTimeIntervalInDay)

    return ES, VaR
