import numpy as np


def HSMeasurements(returns, alpha, weights, portfolioValue, riskMeasureTimeIntervalInDay):
    """
    The function calculates the Value at Risk (VaR) and Expected Shortfall (ES) for a portfolio using the Historical Simulation (HS) approach.

    INPUTS:
        returns:                     The returns of the portfolio assets
        alpha:                       The significance level for the VaR and ES
        weights:                     The weights of the assets in the portfolio
        portfolioValue:              The current value of the portfolio
        riskMeasureTimeIntervalInDay:The time interval for the risk measure in days

    OUTPUTS:
        VaR:                         The Value at Risk for the portfolio
        ES:                          The Expected Shortfall for the portfolio
    """
    nDays = len(returns)
    index = int(nDays * (1 - alpha))
    lossPortfolio = portfolioValue * returns.dot(weights) * np.sqrt(riskMeasureTimeIntervalInDay)  # compute losses of past days
    lossPortfolio = lossPortfolio.sort_values(ascending=False)  # order losses
    VaR = lossPortfolio.iloc[index-1]  # take the m-th largest loss
    ES = np.mean(lossPortfolio.iloc[:index-1])  # compute the mean of losses greater than VaR
    return VaR, ES
