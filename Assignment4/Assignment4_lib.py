
import numpy as np
import pandas as pd
from scipy.stats import t, norm

def AnalyticalNormalMeasures(alpha, nu, weights, portfolioValue, riskMeasureTimeIntervalInDay, returns):

    """
        This function computes the t-student VaR and ES for a given portfolio of assets.
        
        Args:
        - alpha: the confidence level (scalar)
        - nu: the degrees of freedom (scalar)
        - weights: a numpy array of portfolio weights (m, )
        - portfolioValue: the value of the whole portfolio (scalar)
        - riskMeasureTimeIntervalInDay: the time interval for risk measure expressed in days (scalar)
        - returns: a pandas dataframe of asset returns (n, m)
    """

    # find the needed quantile from the student t distribution
    t_alpha = t.ppf(alpha, nu)

    # compute the mean and standard deviation of the loss
    L = - portfolioValue * (returns @ weights)
    mu = np.mean(L)
    sigma = np.std(L)

    # compute the VaR and ES
    VaR = mu + sigma * t_alpha

    # compute the ES
    ES_std = (nu + t_alpha**2) / (nu - 1) * (t.pdf(t_alpha, nu) / (1 - alpha))
    ES = mu + sigma * ES_std

    return ES, VaR

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

def PrincCompAnalysis(yearlyCovariance,  yearlyMeanReturns,  weights,  H,  alpha, numberOfPrincipalComponents, portfolioValue):

    """
        This function computes the non-centered PCA analysis for a given portfolio of assets.
        
        Args:
        - yearlyCovariance: the yearly covariance matrix of the asset returns (m, m)
        - yearlyMeanReturns: the yearly mean returns of the assets (m, )
        - weights: the portfolio weights (m, )
        - H: the time horizon expressed a year fraction (scalar)
        - alpha: the confidence level (scalar)
        - numberOfPrincipalComponents: the number of principal components to consider (scalar)
        - portfolioValue: the value of the whole portfolio (scalar)
    """
    # compute eigenvalues and eigenvectors of the covariance matrix
    eigenvalues, eigenvectors = np.linalg.eig(yearlyCovariance)

    # sort the eigenvalues and eigenvectors
    idx = eigenvalues.argsort()[::-1]
    lambdas = eigenvalues[idx]
    gamma = eigenvectors[:, idx]

    # compute the new weights and mean
    w_hat = gamma.T @ weights
    mu_hat = - gamma.T @ yearlyMeanReturns

    # compute the reduced mean and reduced variance
    mu_red = (w_hat[:numberOfPrincipalComponents] @ mu_hat[:numberOfPrincipalComponents]).sum()
    sigma2_red = (w_hat[:numberOfPrincipalComponents]**2 @ lambdas[:numberOfPrincipalComponents]).sum()

    # compute the VaR and ES
    z_alpha = norm.ppf(alpha)

    VaR = portfolioValue * ( H * mu_red + np.sqrt(H) * np.sqrt(sigma2_red) * z_alpha)

    ES = portfolioValue * ( H * mu_red + np.sqrt(H) * np.sqrt(sigma2_red) * norm.pdf(z_alpha) / (1-alpha))

    # compute the percentage of explained variance
    exp_var_perc = lambdas[:numberOfPrincipalComponents].sum() / lambdas.sum()

    return ES, VaR, exp_var_perc

def FullMonteCarloVaR(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend, \
    volatility, timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha, NumberOfDaysPerYears):

    """
        This function computes the full Monte Carlo VaR for a given portfolio of assets.
        
        Args:
        - logReturns: the log returns of the asset (n, )
        - numberOfShares: the number of shares of the stock (scalar)
        - numberOfCalls: the number of call options (scalar)
        - stockPrice: the current stock price (scalar)
        - strike: the strike price of the put option (scalar)
        - rate: the risk-free rate (scalar)
        - dividend: the dividend yield (scalar)
        - volatility: the volatility of the stock (scalar)
        - timeToMaturityInYears: the time to maturity of the put option expressed in years (scalar)
        - riskMeasureTimeIntervalInYears: the time interval for risk measure expressed in years (scalar)
        - alpha: the confidence level (scalar)
        - NumberOfDaysPerYears: the number of days per year (scalar)
    """

    # define the B&S formula for a call
    def blackScholesCall(S, K, r, d, sigma, T):
        """
    This function computes the price of a call option using the Black-Scholes formula.

    Args:
    - S: the stock price (scalar or numpy array)
    - K: the strike price (scalar)
    - r: the risk-free rate (scalar)
    - d: the dividend yield (scalar)
    - sigma: the volatility of the stock (scalar)
    - T: the time to maturity of the option in years (scalar)
        """
        d1 = (np.log(S/K) + (r - d + 0.5*sigma**2)*T) / (sigma * np.sqrt(T))
        d2 = d1 - sigma * np.sqrt(T)
        return S * np.exp(-d*T) * norm.cdf(d1) - K * np.exp(-r*T) * norm.cdf(d2)

    # compute the price of the call option as of the valuation date
    callPrice = blackScholesCall(stockPrice, strike, rate, dividend, volatility, timeToMaturityInYears)

    # compute the price of the call option as of the valuation date plus the time horizon 
    S_delta = stockPrice * np.exp(logReturns * riskMeasureTimeIntervalInYears * NumberOfDaysPerYears)
    callPrice_delta = blackScholesCall(S_delta, strike, rate, dividend, volatility, timeToMaturityInYears - riskMeasureTimeIntervalInYears)

    # compute the loss as the loss of the stock minus the loss of the call
    loss_call = - numberOfCalls * (callPrice_delta - callPrice)
    loss_stock = - numberOfShares * (S_delta - stockPrice)
    total_loss = loss_stock - loss_call
    
    # take the WHS approach to compute the VaR
    C = (1-alpha)/(1-alpha**len(total_loss))
    weights = [C * alpha**i for i in reversed(range(len(total_loss)))]

    # order the losses and weights and find the VaR
    df_losses = pd.DataFrame({'losses': total_loss, 'weights': weights})
    df_losses = df_losses.sort_values('losses', ascending=False)
    i_star = df_losses[df_losses['weights'].cumsum() <= 1-alpha].index[-1]

    VaR = df_losses.loc[i_star, 'losses']

    return VaR

