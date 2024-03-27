import pandas as pd
import numpy as np
from scipy.stats import t
from scipy.stats import norm
import math

def SliceDataFromStartDate(data, endDate, duration):
    """
    SliceDataFromStartDate takes a dataframe and cuts it above start date of computation of risk measures
    and after the end date
    :param data:          a given data-set
    :param endDate:       the last date wanted
    :param duration:      time window
    :return: data:        the cut data-set
    """

    # reduce dataset up until the end date
    data = data[(data['Date'] <= endDate)]

    # find estimation first date
    endDate = pd.to_datetime(endDate)  # convert it to datetime format
    startDate = endDate - pd.Timedelta(days=duration)  # to go back 5 years in time
    startDate = startDate.strftime('%Y-%m-%d')

    # reduced dataset at the estimation interval
    data = data[(data['Date'] >= startDate)]

    # fill the NaN values with previous values
    data = data.ffill()

    data = data.reset_index(drop=True)

    return data


def AnalyticalNormalMeasures(alpha, weights, portfolioValue, riskMeasureTimeIntervalInDay, returns):
    """
    This function determines VaR and ES of a portfolio using the parametric method
    mean/covariance via t-student with nu=4.

    :param alpha:                           significance value
    :param weights:                         portfolio weights
    :param portfolioValue:                  notional
    :param riskMeasureTimeIntervalInDay:    estimation interval in day
    :param returns:                         returns' matrix

    :return: ES:                            Expected Shortfall
    :return: VaR:                           Value at Risk
    """

    # log returns' mean for each analyzed company
    mu_vector = returns.iloc[:, 1:].mean()
    cov_matrix = pd.DataFrame.cov(returns.iloc[:, 1:])

    # assume loss distribution is gaussian
    mean_loss = - weights.dot(mu_vector)
    cov_loss = (weights.dot(cov_matrix)).dot(weights)

    # degrees of freedom
    nu = 4

    # standard VaR
    VaR_std = t.ppf(alpha, df=nu)

    # VaR with variance/covariance method
    VaR = riskMeasureTimeIntervalInDay * mean_loss + math.sqrt(riskMeasureTimeIntervalInDay) * math.sqrt(
        cov_loss) * VaR_std

    # portfolio VaR
    var_value = portfolioValue * VaR
    print('VaR:', var_value)

    # t-student probability density function
    phi_t = t.pdf(VaR_std, df=nu)

    # standard ES
    ES_std = ((nu + VaR_std ** 2) / (nu - 1)) * (phi_t / (1 - alpha))

    # ES with variance/covariance method
    ES = riskMeasureTimeIntervalInDay * mean_loss + math.sqrt(riskMeasureTimeIntervalInDay) * math.sqrt(
        cov_loss) * ES_std

    # portfolio ES
    es_value = portfolioValue * ES
    print('ES:', es_value)

    return es_value, var_value


def price_to_return(Dataset):
    """
    This function computes log returns from share prices.
    :param Dataset:                           share prices
    :return returns:                          log returns
    """

    # compute the log_returns
    Dataset.iloc[:, 1:] = Dataset.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')
    log_returns = np.log(Dataset.iloc[:, 1:] / Dataset.iloc[:, 1:].shift(1))
    Dataset.iloc[:, 1:] = log_returns
    returns = Dataset.drop(Dataset.index[0])
    returns = returns.reset_index(drop=True)
    return returns


def HSMeasurements(returns, alpha, weights, portfolioValue, riskMeasureTimeIntervalInDay):
    """
    This function performs Historical Simulation measurements of Value at Risk and Expected Shortfall

    :param returns:                         returns' matrix
    :param alpha:                           significance value
    :param weights:                         portfolio weights
    :param portfolioValue:                  notional
    :param riskMeasureTimeIntervalInDay:    estimation interval in day

    :return: ES:                            Expected Shortfall
    :return: VaR:                           Value at Risk
    """

    # "Simulate" loss distribution
    loss = - portfolioValue * returns.iloc[:, 1:].to_numpy().dot(weights)
    # simulated loss distribution using historical log-returns values
    loss.sort(kind='stable')  # sort loss distribution in increasing order
    loss = loss[::-1]  # loss distribution changed to decreasing order (worst loss first value)

    n = len(loss)  # number of historical returns used --> "size" of loss distribution
    index = math.floor(n * (1 - alpha))  # index of the loss corresponding to VaR

    VaR = math.sqrt(riskMeasureTimeIntervalInDay) * loss[index-1]
    print('VaR:', VaR)
    ES = math.sqrt(riskMeasureTimeIntervalInDay) * loss[0:index].mean()
    # mean value of the worst losses up to the VaR one
    print('ES:', ES)

    return ES, VaR


def bootstrapStatistical(numberOfSamplesToBootstrap, returns):
    """
    This function samples M numbers from 1 to n.

    :param numberOfSamplesToBootstrap:      number of samples = M
    :param returns:                         time series of returns of assets in the portfolio
    :return: samples:                       indexes of returns that were sampled
    """
    n = len(returns)  # number of time observations
    samples = np.array([np.random.randint(0, n - 1) for _ in range(numberOfSamplesToBootstrap)])
    return samples


def WHSMeasurements(returns, alpha, lambda_P, weights, portfolioValue, riskMeasureTimeIntervalInDay):
    """
    This function determines VaR and ES of a portfolio using Weighted Historical Simulation

    :param returns:                         returns' matrix
    :param alpha:                           significance value
    :param lambda_P:                        historical index
    :param weights:                         portfolio weights
    :param portfolioValue:                  notional
    :param riskMeasureTimeIntervalInDay:    estimation interval in day
    :return: ES:                            Expected Shortfall with WHS
    :return: VaR:                           Value at Risk with WHS
    """

    # observations number
    n = len(returns)

    # normalization factor
    C = (1 - lambda_P) / (1 - lambda_P ** n)

    # historical losses
    L = - portfolioValue * returns.iloc[:, 1:].to_numpy().dot(weights)

    # Compute the exponents for the lambda coefficients: decreasing in time
    lambdaExponent = np.arange(n - 1, -1, -1)

    # compute simulation weights
    weights_sim = C * np.power(lambda_P, lambdaExponent)

    # order losses in decreasing way and the respective simulation weights
    L, weights_sim = zip(*sorted(zip(L, weights_sim), reverse=True))
    # EXPLANATION OF ZIP(*SORTED(ZIP...)): the first zip binds together L and weights_sim as sort of pairs of
    # coordinates (x,y); then *sorted orders these pairs based on the ordering of the x value (so the losses L), and we
    # specified reverse=True to have the losses ordered from biggest to smallest; finally, the second zip "unzips" the
    # pairs and returns them as two separate vectors, but now they're ordered how we wanted and the relationship
    # between losses and historical weights are maintained

    # find the index that satisfy the constraints
    # initialize index counter
    i_temp = 0
    # initialize sum simulation weights
    sum_weights_sim = 0

    # while loop to find the index of the weights sum corresponding to 1-alpha
    while sum_weights_sim <= (1 - alpha):
        sum_weights_sim += weights_sim[i_temp]
        i_temp += 1
    i_star = i_temp - 1

    # compute VaR with WHS
    VaR_WHS = math.sqrt(riskMeasureTimeIntervalInDay) * L[i_star]
    print('VaR:', VaR_WHS)

    # compute ES with WHS
    ES_WHS = (math.sqrt(riskMeasureTimeIntervalInDay) * (np.dot(weights_sim[:i_star + 1], L[:i_star + 1])) /
              (np.sum(weights_sim[:i_star + 1])))
    print('ES:', ES_WHS)

    return ES_WHS, VaR_WHS


def plausibilityCheck(returns, portfolioWeights, alpha, portfolioValue, riskMeasureTimeIntervalInDay):
    """
    This function performs a plausibility check on the VaR (rule of thumb)

    :param returns:                         returns' matrix
    :param portfolioWeights:                weights of portfolio = sensitivities
    :param alpha:                           significance value
    :param portfolioValue:                  notional
    :param riskMeasureTimeIntervalInDay:    estimation interval in day
    :return: VaR:                           Value at Risk
    """
    # Remove the date from the returns dataFrame and transform it to a numpy array
    returns = returns.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')
    returns = returns.reset_index(drop=True)

    # Compute upper and lower quantiles of the risk factors distributions (asset returns), one for each asset
    upperQuantiles = np.quantile(returns, alpha, axis=0)  # axis = 0 to compute one quantile for each column
    lowerQuantiles = np.quantile(returns, 1 - alpha, axis=0)

    # Compute signedVaR for each risk factor
    signedVaR = portfolioWeights * (abs(upperQuantiles) + abs(lowerQuantiles)) / 2

    # Compute Correlation Matrix with command corrcoef
    CorrelationMatrix = np.corrcoef(returns, rowvar=False)  # rowvar = False means that each column is a variable

    # Compute VaR according to plausibility check rule of thumb
    VaR = (math.sqrt(signedVaR.transpose().dot(CorrelationMatrix).dot(signedVaR)) *
           math.sqrt(riskMeasureTimeIntervalInDay))
    VaR_value = VaR * portfolioValue  # multiply for the portfolio value in time t to get reasonable measure
    print('Plausibility check on VaR:', VaR_value)

    return VaR_value


def PrincCompAnalysis(yearlyCovariance, yearlyMeanReturns, weights, H, alpha, numberOfPrincipalComponents,
                      portfolioValue):
    """
    this function calculates the VaR and the ES via a Gaussian parametric PCA approach

    :param yearlyCovariance:                matrix of yearly covariances of returns
    :param yearlyMeanReturns:               vector of yearly means of returns
    :param weights:                         weights of portfolio
    :param H:                               estimation interval in fraction of years
    :param alpha:                           confidence level
    :param numberOfPrincipalComponents:     number of principal components to be analyzed
    :param portfolioValue:                  value of portfolio
    :return: ES:                            Expected Shortfall
    :return: VaR:                           Value at Risk
    """
    # Compute VaR std and ES std for Gaussian distributed losses (as by assumption of PCA method)
    VaR_std = norm.ppf(alpha, loc=0, scale=1)
    ES_std = norm.cdf(norm.ppf(alpha, loc=0, scale=1)) / (1 - alpha)

    # compute eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eig(yearlyCovariance)
    # order for decreasing eigenvalues
    eigenvalues, eigenvectors = zip(*sorted(zip(eigenvalues, eigenvectors.T), reverse=True))
    # zip functioning explained in WHS method
    eigenvalues = np.diag(eigenvalues)  # make eig into a diag matrix
    eigenvectors = np.array(eigenvectors).T  # make eigenvectors be  matrix again, it was messed up cause of zip

    # new definitions
    mu_hat = np.dot(eigenvectors.T, yearlyMeanReturns)
    weights_hat = np.dot(eigenvectors.T, weights)

    # mean and covariance for reduced form portfolio
    mu_rfp = sum((weights_hat * mu_hat)[:numberOfPrincipalComponents])
    sigma_rfp = sum(((weights_hat ** 2).dot(eigenvalues))[:numberOfPrincipalComponents])

    # computation of VaR
    VaR = H * mu_rfp + math.sqrt(H) * math.sqrt(sigma_rfp) * VaR_std
    VaR_value = portfolioValue * VaR

    # computation of ES
    ES = H * mu_rfp + math.sqrt(H) * math.sqrt(sigma_rfp) * ES_std
    ES_value = ES * portfolioValue

    return ES_value, VaR_value


def FullMonteCarloVaR(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend, volatility,
                      timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha, NumberOfDaysPerYears, lambdaWHS):
    """
    Function to compute full Monte Carlo VaR of portfolio

    :param logReturns:                          matrix of log-returns
    :param numberOfShares:                      number of shares of stock held in portfolio
    :param numberOfCalls:                       number of calls held in portfolio
    :param stockPrice:                          price of stock
    :param strike:                              strike value of calls
    :param rate:                                fixed interest rate of stock
    :param dividend:                            dividend yield of stock
    :param volatility:                          volatility of stock
    :param timeToMaturityInYears:               time to maturity of calls
    :param riskMeasureTimeIntervalInYears:      estimation interval expressed in fraction of years
    :param alpha:                               confidence level
    :param NumberOfDaysPerYears:                business days considered per year
    :param lambdaWHS:                           parameter for Weighted Historical Simulation
    :return:VaR_WHS:                            Value at Risk
    """
    # random.seed(40)  # for reproducibility of results reasons

    # Black-Scholes formula for Call price
    callPrice = blsCall(stockPrice, strike, rate, dividend, volatility, timeToMaturityInYears)

    # historical log returns in the past over the 10 days time lag
    logReturns = logReturns.iloc[:, 1].to_numpy()
    delta = int(riskMeasureTimeIntervalInYears * NumberOfDaysPerYears)
    nObs = int(len(logReturns))

    # Randomly extract indexes from the vector of logReturns to use as sampled log returns
    nSim = int(1e5)

    # Parameter C for the weights of weighted historical simulation
    C = (1 - lambdaWHS) / (1 - lambdaWHS ** nObs)

    # Compute the weight for each index of our past log-returns, using the WHS rule
    initialIdxVec = np.arange(0, nObs, 1)
    initialWeights = C * lambdaWHS ** np.arange(nObs - 1, -1, -1)

    # Extract 10 indexes for each simulation, not extracting from uniform but giving more prob. to the higher weights
    indexes = np.random.choice(initialIdxVec, size=(nSim, 10), p=initialWeights)
    # we want 10-days VaR, so we need 10 days log-returns, so for each simulation I take 10 indexes, so I simulate
    # 10-day returns summing 10 random daily returns

    # Matrix of log returns as taken above
    logRetMC = logReturns[indexes]  # extract sampled log returns
    logRetDeltaMC = np.sum(logRetMC, axis=1)  # sum the values of the ten log-returns for each simulation

    # Find the weights of each simulation by taking the mean index for each simulation and computing its relative weight
    meanIndexes = np.around(np.mean(indexes, axis=1)).astype('i')  # mean for each row = or each simulation
    weightsSims = initialWeights[meanIndexes]  # extract relative weights for each index
    weightsSims = weightsSims/np.sum(weightsSims)  # normalize the weights so they sum up to 1

    # Evaluate vector of new prices in t+delta: one for each simulation
    vectorStockt1 = stockPrice * np.exp(logRetDeltaMC)

    # Evaluate call prices considering the simulated vector of prices in t+delta
    timeToMaturityInYears -= delta/365
    vectorCallt1 = blsCall(vectorStockt1, strike, rate, dividend, volatility, timeToMaturityInYears)

    # Evaluate derivative losses (Remark: we're shorting the derivative, so we remove the minus)
    lossDer = numberOfCalls * (vectorCallt1 - callPrice)

    # Evaluate stock losses
    lossStock = - numberOfShares * (vectorStockt1 - stockPrice)

    # Loss total portfolio
    lossTotal = np.add(lossDer, lossStock)

    # Now I have vector of total losses for each simulation, and I also have the corresponding weights of each
    # simulation stored in the weightsSims vector
    # Sorting the losses in decreasing order, using zip sorted (keeping relation with weights)
    lossTotal, weightsSims = zip(*sorted(zip(lossTotal, weightsSims), reverse=True))
    # zip functioning explained in WHS method

    # find the index that satisfy the constraints
    # initialize index counter
    i_temp = 0
    # initialize sum simulation weights
    sum_weights_sim = 0
    # while loop to find the index of the weights sum corresponding to 1-alpha
    while sum_weights_sim <= (1 - alpha):
        sum_weights_sim += weightsSims[i_temp]
        i_temp += 1
    i_star = i_temp - 1

    # compute VaR with WHS
    VaR = lossTotal[i_star]
    print('VaR:', float(VaR))

    return VaR


def DeltaNormalVaR(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend, volatility,
                   timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha, NumberOfDaysPerYears, lambdaWHS):
    """
    Function to compute full VaR of portfolio with Delta Normal method

    :param logReturns:                          matrix of log-returns
    :param numberOfShares:                      number of shares of stock held in portfolio
    :param numberOfCalls:                       number of calls held in portfolio
    :param stockPrice:                          price of stock
    :param strike:                              strike value of calls
    :param rate:                                fixed interest rate of stock
    :param dividend:                            dividend yield of stock
    :param volatility:                          volatility of stock
    :param timeToMaturityInYears:               time to maturity of calls
    :param riskMeasureTimeIntervalInYears:      estimation interval expressed in fraction of years
    :param alpha:                               confidence level
    :param NumberOfDaysPerYears:                business days considered per year
    :param lambdaWHS:                           parameter for Weighted Historical Simulation
    :return:VaR_WHS:                            Value at Risk
    """

    # historical log returns in the past over the 10 days time lag
    logReturns = logReturns.iloc[:, 1].to_numpy()
    # num observations
    nObs = len(logReturns)
    delta = int(riskMeasureTimeIntervalInYears * NumberOfDaysPerYears)
    nSim = int(1e4)

    # Parameter C for the weights of weighted historical simulation
    C = (1 - lambdaWHS) / (1 - lambdaWHS ** nObs)

    # Compute the weight for each index of our past log-returns, using the WHS rule
    initialIdxVec = np.arange(0, nObs, 1)
    initialWeights = C * lambdaWHS ** np.arange(nObs - 1, -1, -1)

    # Extract 10 indexes for each simulation, not extracting from uniform but giving more prob. to the higher weights
    indexes = np.random.choice(initialIdxVec, size=(nSim, 10), p=initialWeights)
    # we want 10-days VaR, so we need 10 days log-returns, so for each simulation I take 10 indexes, so I simulate
    # 10-day returns summing 10 random daily returns

    # Matrix of log returns as taken above
    logRetMC = logReturns[indexes]  # extract sampled log returns
    logRetDeltaMC = np.sum(logRetMC, axis=1)  # sum the values of the ten log-returns for each simulation

    # Find the weights of each simulation by taking the mean index for each simulation and computing its relative weight
    meanIndexes = np.around(np.mean(indexes, axis=1)).astype('i')  # mean for each row = or each simulation
    weightsSims = initialWeights[meanIndexes]  # extract relative weights for each index
    weightsSims = weightsSims/np.sum(weightsSims)  # normalize the weights so they sum up to 1

    # Evaluate sensitivity of the Call
    d1 = (np.log(stockPrice / strike) + (rate - dividend + 0.5 * volatility ** 2) * timeToMaturityInYears) / (
            volatility * np.sqrt(timeToMaturityInYears))
    # Delta for Call with given params
    deltaCall = np.exp(-dividend * timeToMaturityInYears) * norm.cdf(d1)

    # Loss total portfolio
    lossTotal = -(-numberOfCalls * deltaCall + numberOfShares) * stockPrice * logRetDeltaMC
    # Now I have a vector of total losses, which simulates the loss distribution, and I also have
    # the corresponding weights of each observation stored in the weightsSims vector

    # Sorting the losses in decreasing order, using zip sorted (keeping relation with weights)
    lossTotal, weightsSims = zip(*sorted(zip(lossTotal, weightsSims), reverse=True))
    # zip functioning explained in WHS method

    # find the index that satisfy the constraints
    # initialize index counter
    i_temp = 0
    # initialize sum simulation weights
    sum_weights_sim = 0

    # while loop to find the index of the weights sum corresponding to 1-alpha
    while sum_weights_sim <= (1 - alpha):
        sum_weights_sim += weightsSims[i_temp]
        i_temp += 1
    i_star = i_temp - 1

    # compute VaR with WHS, using scaling factor because we computed a daily VaR and we want to have a 10day VaR
    VaR = lossTotal[i_star]
    print('VaR:', float(VaR))

    return VaR


def blsCall(stockPrice, strike, rate, dividend, volatility, timeToMaturityInYears):
    """
    Calculates the Call option price according to Black and Scholes formula

    :param stockPrice:              today's price of stock
    :param strike:                  strike value of call
    :param rate:                    fixed interest rate
    :param dividend:                dividend yield of stock
    :param volatility:              volatility of stock
    :param timeToMaturityInYears:   time to maturity of call
    :return: callPrice:             price of call option
    """

    # Black-Scholes parameters d1 and d2
    d1 = (np.log(stockPrice / strike) + (rate - dividend + 0.5 * volatility ** 2) * timeToMaturityInYears) / (
            volatility * np.sqrt(timeToMaturityInYears))
    d2 = d1 - volatility * np.sqrt(timeToMaturityInYears)

    # Evaluate call prices
    callPrice = stockPrice * np.exp(-dividend * timeToMaturityInYears) * norm.cdf(d1) - strike * np.exp(
        -rate * timeToMaturityInYears) * norm.cdf(d2)

    return callPrice
