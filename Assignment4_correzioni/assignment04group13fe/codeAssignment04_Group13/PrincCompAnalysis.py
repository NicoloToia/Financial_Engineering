import numpy as np
from scipy.stats import norm


def PrincCompAnalysis(yearlyCovariance, yearlyMeanReturns, weights, H, alpha, numberOfPrincipalComponents, portfolioValue):
    """
    The function performs Principal Component Analysis (PCA) on a portfolio and calculates the Value at Risk (VaR) and Expected Shortfall (ES).

    INPUTS:
        yearlyCovariance:             The yearly covariance of the portfolio returns
        yearlyMeanReturns:            The yearly mean returns of the portfolio
        weights:                      The weights of the assets in the portfolio
        H:                            The horizon for the risk measure
        alpha:                        The significance level for the VaR and ES
        numberOfPrincipalComponents:  The number of principal components to consider in the PCA
        portfolioValue:               The current value of the portfolio

    OUTPUTS:
        VaR:  The Value at Risk for the portfolio
        ES:   The Expected Shortfall for the portfolio
    """
    n = numberOfPrincipalComponents
    mu = - yearlyMeanReturns
    Sigma = yearlyCovariance
    eigenvalues, eigenvectors = np.linalg.eig(Sigma)  # eigenvalues and eigenvectors of covariance matrix
    sorted_indices = np.argsort(eigenvalues)[::-1]  # obtain indices of ordered eigenvalues
    sorted_eigenvalues = eigenvalues[sorted_indices]  # order eigenvalues
    sorted_eigenvectors = eigenvectors[:, sorted_indices]  # order eigenvectors based on ordered eigenvalues
    Gamma = sorted_eigenvectors
    mu_sorted = mu.iloc[sorted_indices]  # order means based on ordered eigenvalues
    mu_hat = np.dot(Gamma.T, mu_sorted)
    omega_hat = np.dot(Gamma.T, weights)
    mu_red = mu_hat[:n].dot(omega_hat[:n])  # compute reduced means
    var_red = (omega_hat[:n]**2).dot(sorted_eigenvalues[:n])  # compute reduced variances
    VaR = H * mu_red + np.sqrt(H) * np.sqrt(var_red) * norm.ppf(alpha)  # compute VaR of the PCA
    ES = H * mu_red + np.sqrt(H) * np.sqrt(var_red) * (norm.pdf(norm.ppf(alpha))/(1 - alpha))  # compute ES of the PCA
    return VaR * portfolioValue, ES * portfolioValue
