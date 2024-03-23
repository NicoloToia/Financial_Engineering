
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

def WHSMeasurements(returns,  alpha,  lambda_,  weights,  portfolioValue, riskMeasureTimeIntervalInDay):

    """
        This function computes the weighted historical simulation VaR and ES for a given portfolio of assets.

        Args:
        - returns: a pandas dataframe of asset returns (n, m)
        - alpha: the confidence level (scalar)
        - lambda_ : the decay factor (scalar) (cannot use lambda as it is a reserved keyword in Python)
        - weights: a numpy array of portfolio weights (m, )
        - portfolioValue: the value of the whole portfolio (scalar)
        - riskMeasureTimeIntervalInDay: the time interval for risk measure expressed in days (scalar)

    """
    
    losses = - portfolioValue * (returns @ weights)

    # compute the corresponding time weights
    C = (1-lambda_)/(1-lambda_**len(losses)) # normalization factor
    weights = [C * lambda_**i for i in reversed(range(len(losses)))]

    # make the losses and weights into a dataframe and sort it by losses (greatest first)
    df_losses = pd.DataFrame({'losses': losses, 'weights': weights})
    df_losses = df_losses.sort_values('losses', ascending=False)

    # find the index where the cumulative sum of the weights is lower or equal to 1-alpha
    i_star = df_losses[df_losses['weights'].cumsum() <= 1-alpha].index[-1]

    # VaR is the corresponding loss
    VaR = np.sqrt(riskMeasureTimeIntervalInDay) * df_losses.loc[i_star, 'losses']
    # ES is the weighted average of the losses
    sum_weights = df_losses.loc[:i_star, 'weights'].sum()
    sum_losses_weighted = (df_losses.loc[:i_star, 'losses'] * df_losses.loc[:i_star, 'weights']).sum()

    ES = np.sqrt(riskMeasureTimeIntervalInDay) * sum_losses_weighted / sum_weights

    return ES, VaR

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
    lb = returns.quantile(1-alpha, interpolation='nearest')
    ub = returns.quantile(alpha, interpolation='nearest')

    # compute the signed VaR
    sVaR = portfolioValue * portfolioWeights * (abs(lb) + abs(ub)) / 2

    # correlation matrix
    corr = returns.corr()

    # compute the VaR with the thumb rule
    VaR = np.sqrt(riskMeasureTimeIntervalInDay) * np.sqrt(sVaR @ corr @ sVaR)

    return VaR