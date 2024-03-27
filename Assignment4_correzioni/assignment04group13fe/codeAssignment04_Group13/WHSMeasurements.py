import numpy as np


def WHSMeasurements(returns, alpha, lambda_, weights, portfolioValue, riskMeasureTimeIntervalInDay):
    """
    The function calculates the Value at Risk (VaR) and Expected Shortfall (ES) for a portfolio using the Weighted Historical Simulation (WHS) approach.

    INPUTS:
        returns:                     The returns of the portfolio assets
        alpha:                       The significance level for the VaR and ES
        lambda_:                     The decay factor for the weights in the WHS
        weights:                     The weights of the assets in the portfolio
        portfolioValue:              The current value of the portfolio
        riskMeasureTimeIntervalInDay:The time interval for the risk measure in days

    OUTPUTS:
        VaR:                         The Value at Risk for the portfolio
        ES:                          The Expected Shortfall for the portfolio
    """
    C = (1 - lambda_)/(1 - lambda_ ** len(returns))
    s = np.arange(1, len(returns)+1, 1)
    t = len(returns)
    omega = C * lambda_ ** (t - s)  # weights of past losses
    lossPortfolio = - portfolioValue * returns.dot(weights) * np.sqrt(riskMeasureTimeIntervalInDay)  # compute losses of past days
    lossSortedIndices = np.argsort(- lossPortfolio)  # obtain indices of ordered losses
    sortedOmega = omega[lossSortedIndices]  # order weights based on ordered losses
    sortedLoss = lossPortfolio.sort_values(ascending=False)  # order losses
    iStar = searchIndex(sortedOmega, alpha)
    VaR = sortedLoss.iloc[iStar]  # take the iStar-th largest loss
    ES = (sortedOmega[:iStar].dot(sortedLoss.iloc[:iStar])/np.sum(sortedOmega[:iStar]))  # compute the mean of losses greater than VaR
    return VaR, ES


# this function finds the minimum index for that the condition of the WHSM is satisfied
def searchIndex(omega, alpha):
    """
    The function finds the minimum index for which the condition of the WHSM is satisfied.

    INPUTS:
        omega:  The weights of past losses
        alpha:  The significance level for the VaR and ES

    OUTPUTS:
        iStar:  The minimum index for which the condition of the WHSM is satisfied
    """
    sum = np.cumsum(omega)
    iStar = np.searchsorted(sum, 1 - alpha, side='right')
    return iStar

