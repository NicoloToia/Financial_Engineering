import numpy as np


def plausibilityCheck(returns, portfolioWeights, alpha, portfolioValue, riskMeasureTimeIntervalInDay):
    """
    The function performs a plausibility check on the portfolio by calculating the Value at Risk (VaR) using the correlation coefficients and percentiles of the returns.

    INPUTS:
        returns:                     The returns of the portfolio assets
        portfolioWeights:            The weights of the assets in the portfolio
        alpha:                       The significance level for the VaR
        portfolioValue:              The current value of the portfolio
        riskMeasureTimeIntervalInDay:The time interval for the risk measure in days

    OUTPUTS:
        VaR:                         The Value at Risk for the portfolio
    """
    C = np.corrcoef(returns, rowvar=False)  # correlation coefficients
    l = np.percentile(returns, (1 - alpha) * 100, axis=0)  # lower percentile
    u = np.percentile(returns, alpha * 100, axis=0)  # upper percentile
    sens = - portfolioValue * portfolioWeights  # compute sensitivities (linear portfolio)
    sVaR = sens * (abs(l) + abs(u)) / 2
    VaR = np.sqrt(riskMeasureTimeIntervalInDay) * np.sqrt(np.dot(sVaR.T, np.dot(C, sVaR)))
    return VaR
