
# functions we need
import numpy as np
import pandas as pd
import scipy as sc
import datetime
import math
import random

def plausibilityCheck(returns, portfolioWeights, alpha, portfolioValue, riskMeasureTimeIntervalInDay):

  ## The function compute the plausability check in the VaR computation for a stock portfolio##

  # REMARK ON INPUT : returns is the dataframe of log returns

  # Taking just the values (without indexes)
  log_returns = returns.values
  # Computing the portfolio sensitivity
  sens = -portfolioValue*portfolioWeights;

  # Computing the quantile of the risk factors
  l_down = np.quantile(log_returns,1-alpha,0)
  l_up = np.quantile(log_returns,alpha,0)

  # Computing the signed-VaR
  mean_value = (np.abs(l_down)+np.abs(l_up))/2
  sVaR = sens*mean_value

  # Computing the correlation matrix
  sigma = returns.corr()
  sigma = sigma.to_numpy()

  # Computing the VaR
  VaR = np.sqrt(np.dot(np.dot(sVaR,sigma),sVaR.T))*math.sqrt(riskMeasureTimeIntervalInDay)

  return VaR

def AnalyticalNormalMeasures(alpha, weights, portfolioValue, riskMeasureTimeIntervalInDay, returns):

  ## Function for computing the VaR with variance-covariance method with t_student (nu=4) for a stock portfolio##

  # REMARK ON INPUT : returns are log returns

  # Computing the mean as numpy array
  mu = returns.mean()
  mu = mu.to_numpy()

  # Variance-Covariance matrix
  sigma = returns.cov()
  sigma = sigma.to_numpy()

  # Computing the loss mean and variance for the portfolio with unitary value
  mu_loss = - np.dot(weights,mu)
  sigma_loss = np.dot(np.dot(weights,sigma),weights)

  # Computing the quantile for the t-student with 4 degrees of freedom
  nu = 4
  quant = sc.stats.t.ppf(alpha,nu)

  # Computing the VaR
  VaR = riskMeasureTimeIntervalInDay*mu_loss + math.sqrt(riskMeasureTimeIntervalInDay)*math.sqrt(sigma_loss)*quant
  VaR = VaR*portfolioValue

  # Computing the standardized ES for a t-student
  ES_std = (nu + quant**2)/(nu-1) * (sc.stats.t.pdf(quant,nu))/(1-alpha)

  # Computing the ES
  ES = (riskMeasureTimeIntervalInDay*mu_loss + math.sqrt(riskMeasureTimeIntervalInDay)*math.sqrt(sigma_loss)*ES_std)*portfolioValue

  return [VaR,ES]

def HSMeasurements(returns, alpha, weights, portfolioValue, riskMeasureTimeIntervalInDay):

  ## Function for computing VaR via Historical Simulation method for a stock portfolio##

  # REMARK ON INPUT : returns is the dataframe of log returns

  # Taking log returns as numpy matrix
  log_ret = returns.values

  # Computing the loss for each past date
  loss = - portfolioValue*np.dot(log_ret,weights)

  # Sort the loss in decreasing order
  loss.sort()
  loss = np.flip(loss)

  # Computing the VaR as the loss vector quantile
  business_days = returns.shape[0]-1
  VaR = loss[int(business_days*(1-alpha))+1]*math.sqrt(riskMeasureTimeIntervalInDay)

  # Computing the ES as the mean of the losses greater than the VaR
  ES = np.mean(loss[1:int(business_days*(1-alpha)+2)])*math.sqrt(riskMeasureTimeIntervalInDay)

  return [VaR, ES]

def bootstrapStatistical(numberOfSamplesToBootstrap, returns):

  ## Function for sampling a prespecified number (given as input) of log returns ##

  # REMARK ON INPUT : returns must be a dataframe (with dates as index)

  # The sample function avoid repetition
  samples = returns.sample(n=numberOfSamplesToBootstrap, replace=True)

  return samples

def WHSMeasurements(returns, alpha, lamb, weights, portfolioValue, riskMeasureTimeIntervalInDay):

  ## The function implements the Weighted Historical Simulation approach for computing VaR and Expected Shortfall for a stock porfolio ##

  # REMARK ON INPUT : returns is the dataframe of log returns

  # Taking log returns as numpy matrix
  log_ret = returns.values

  # Computing the loss in each past date
  loss = -portfolioValue*np.dot(log_ret,weights)

  # Computing the weights in decreasing order (recent dates must be more relevant)
  n = log_ret.shape[0]-1
  c = (1-lamb)/(1-lamb**n)
  loss_weights = [c*lamb**i for i in range(n)]
  loss_weights = np.flip(loss_weights)

  # Linking every loss with its weight and sorting them by the loss
  combined = list(zip(loss, loss_weights))
  sorted_combined = sorted(combined, key=lambda x: x[0]) # Increasing losses

  # Ordering losses and weights in decreasing order
  loss = np.flip([item[0] for item in sorted_combined])
  loss_weights = np.flip([item[1] for item in sorted_combined])

  # Computing the quantile as the largest loss such that the cumulative weight of the loss is smaller or equal than 1 - alpha
  sum = 0
  cont = -1
  while sum<=1-alpha:
    cont+=1
    sum+=loss_weights[cont]
  cont-=1

  # Computing VaR and ES
  VaR = loss[cont]*math.sqrt(riskMeasureTimeIntervalInDay)
  ES = (np.dot(loss_weights[0:cont+1],loss[0:cont+1])/np.sum(loss_weights[0:cont+1]))*math.sqrt(riskMeasureTimeIntervalInDay)

  return [VaR, ES]

def PrincCompAnalysis(yearlyCovariance, yearlyMeanReturns, weights, H, alpha, numberOfPrincipalComponents, portfolioValue):

  ## The function evaluates VaR  for a stock portfolio via Principal Components Analysis ##

  # Computing the eigenvalues and eigenvector of the covariance matrix with the function linalg.eig of numpy
  eigenval,eigenvect = np.linalg.eig(yearlyCovariance)

  # Computing the indexes of the eigenvalues in decreasing order
  decreasing_indexes = np.argsort(eigenval)[::-1]

  # Computing the gamma matrix with the eigenvectors in the right order
  gamma = eigenvect[:, decreasing_indexes]

  # Reordering eigenvalues in decreasing order
  decreasing_eigenval = eigenval[decreasing_indexes]
  decreasing_eigenval = decreasing_eigenval.reshape(18,1)

  # Computing the projection for the mean vector and the weights through gamma
  mu_cap = np.dot(gamma.T,yearlyMeanReturns)
  w_cap = np.dot(gamma.T,weights.T)
  w_cap_squared = w_cap * w_cap

  # Computing mean and variance for the reduced portfolio
  sigma_red_squared = np.sum(w_cap_squared[0:numberOfPrincipalComponents]*decreasing_eigenval[0:numberOfPrincipalComponents])
  mu_red = np.sum(w_cap[0:numberOfPrincipalComponents]*mu_cap[0:numberOfPrincipalComponents])

  # Computing the quantile of order alpha for a gaussian random variable
  quant = sc.stats.norm.ppf(alpha)

  # Computing VaR
  VaR = H*mu_red + math.sqrt(H)*np.sqrt(sigma_red_squared)*quant #VaR

  # Computing ES
  ES_std = sc.stats.norm.pdf(quant)/(1-alpha)
  ES = H*mu_red + math.sqrt(H)*np.sqrt(sigma_red_squared)*ES_std #ES

  return [VaR, ES]

def BS_CALL(S, K, T, r, q, sigma):
  # The function computes the European Call Option price
  d1 = (np.log(S/K) + (r - q + (sigma**2)/2)*T) / (sigma*np.sqrt(T))
  d2 = d1 - sigma* np.sqrt(T)
  return S*np.exp(-q*T) * sc.stats.norm.cdf(d1) - K * np.exp(-r*T)* sc.stats.norm.cdf(d2)

def DELTA_CALL(S,K,TTM,r,q,sigma):
  # The function computes the delta for an European Call Option
  d1 = (np.log(S/K) + (r - q + (sigma**2)/2)*TTM) / (sigma*np.sqrt(TTM))
  return sc.stats.norm.cdf(d1)*np.exp(-q*TTM)

def FullMonteCarloVaR2(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend, volatility, timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha,NumberOfDaysPerYears):

  ## The Function computes VaR via full Monte Carlo evaluation for a portfolio composed by stocks of a single firm and calls having the same stock as underlying ##

  # The procedure consists on weighting as in the WSH approach the log returns and on the random sampling of a prespecified number of log returns among them
  # These log returns are used for the simulation of the underlying, and then thanks to the closed formula we can price the derivative
  # We are now able to compute the portfolio loss, and the VaR following a WHS approach
  # These operations (very similar to a statistical bootstrap approach) are repeated N times (in the function the variable name is numIter), obtaining N VaRs
  # The final VaR is computed as the mean of the ones computed

  # REMARK ON INPUT : returns is the dataframe of log returns

  # Computing the WHS weights
  lamb = 0.95
  n = logReturns.shape[0]
  c = (1-lamb)/(1-lamb**n)
  weights = [c*lamb**i for i in range(n)]
  weights = np.flip(weights)

  # Linking logreturns and weights in a dataframe
  logReturns["weights"] = weights

  # Setting the number of iteration
  numIter = 1000

  # Setting the number of sample
  extract = 450
  # In this exercise we use 514 returns, so extract can not be greater than 514

  # initializing the vector for the VaRs
  VaR = np.zeros((numIter,1))

  # Computing the initial portfolio value
  initialPtfValue = stockPrice*numberOfShares-numberOfCalls*BS_CALL(stockPrice,strike,timeToMaturityInYears,rate,dividend,volatility)

  # Populating the VaRs vector
  for i in range(numIter):
    # Sampling random log returns from our data
    sample = logReturns.sample(n=extract)
    # Computing the simulated stock value
    stockValue = stockPrice*np.exp(sample.values[:,0])
    # Pricing the call using the simulated stock value
    price_call = BS_CALL(stockValue,strike,timeToMaturityInYears-1/NumberOfDaysPerYears,rate,dividend,volatility)
    # Computing the portfolio loss
    loss = initialPtfValue - (numberOfShares*stockValue-numberOfCalls*price_call)
    # Linking loss and the dataframe logReturns which contains also the weights
    combined = list(zip(loss, sample.values[:,1]))
    # Sorting the loss and flipping losses and weights in decreasing order
    sorted_combined = sorted(combined, key=lambda x: x[0])
    loss = np.flip([item[0] for item in sorted_combined])
    loss_weights = np.flip([item[1] for item in sorted_combined])
    # Computing the quantile with an approach "WHS style"
    sum = 0
    cont = -1
    while sum<=1-alpha:
      cont+=1
      sum+=loss_weights[cont]
    cont-=1
    # Computing the VaR in the single iteration
    VaR[i] = loss[cont]

  # Evaluating the VaR as mean of the VaR
  VaR_MC = np.sum(VaR)/numIter*math.sqrt(riskMeasureTimeIntervalInYears*NumberOfDaysPerYears)

  VaR_medio = np.mean(VaR)*math.sqrt(riskMeasureTimeIntervalInYears*NumberOfDaysPerYears)
  VaR_std = np.std(VaR)*math.sqrt(riskMeasureTimeIntervalInYears*NumberOfDaysPerYears)
  quant = sc.stats.norm.ppf(alpha)
  return VaR_medio-VaR_std*quant, VaR_medio + VaR_std*quant

def FullMonteCarloVaR(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend, volatility, timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha,NumberOfDaysPerYears):

  ## The function compute the VaR with the full evaluation Monte Carlo method ##

  # REMARK ON INPUT : returns is the dataframe of log returns

  # Computing the weights with WHS approach
  lamb = 0.95
  n = logReturns.shape[0]
  c = (1-lamb)/(1-lamb**n)
  weights = [c*lamb**i for i in range(n)]
  pesi = np.flip(weights)

  # Setting the number of sample
  numIter = 10000

  # Initial portfolio value
  initialPtfValue = stockPrice*numberOfShares-numberOfCalls*BS_CALL(stockPrice,strike,timeToMaturityInYears,rate,dividend,volatility)

  # initializing the loss vector
  loss = np.zeros((numIter,1))

  # MC simulation
  for i in range(numIter):
    samples = random.choices(logReturns.values, weights=pesi, k=10)
    return_10gg = np.sum(samples)
    future_stock = stockPrice*np.exp(return_10gg)
    future_call = BS_CALL(future_stock,strike,timeToMaturityInYears-riskMeasureTimeIntervalInYears/NumberOfDaysPerYears,rate,dividend,volatility)
    loss[i] = initialPtfValue - future_stock*numberOfShares + numberOfCalls*future_call

  # Taking the quantile of the loss distribution
  VaR = np.quantile(loss,alpha)

  return VaR

def DeltaNormalVaR(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend, volatility, timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha, NumberOfDaysPerYears):

  ## The function implements the Delta-Normal method for computing the VaR for a portfolio made of stocks of a specific firm and European Call Options written on that underlying ##

 # REMARK ON INPUT : returns is the dataframe of log returns

  # Computating log returns' weights (WHS approach)
  lamb = 0.95
  n = logReturns.shape[0]
  c = (1-lamb)/(1-lamb**n)
  weights = [c*lamb**i for i in range(n)]
  weights = np.flip(weights)

  future_stock = stockPrice*np.exp(logReturns.values) #future value of the stock ("tomorrow's value")

  # Loss of the stock
  loss_stock = - numberOfShares*stockPrice*logReturns.values

  # Loss of the call
  loss_call = - numberOfCalls*DELTA_CALL(future_stock,strike,timeToMaturityInYears-1/NumberOfDaysPerYears,rate,dividend,volatility)*stockPrice*logReturns.values

  # Loss of the entire ptf, we use "-" since we are shortselling the options
  loss = loss_stock - loss_call

  # Sorting losses and the associated weights
  combined = list(zip(loss, weights))
  sorted_combined = sorted(combined, key=lambda x: x[0])

  # We want losses in decreasing order
  loss = np.flip([item[0] for item in sorted_combined])
  loss_weights = np.flip([item[1] for item in sorted_combined])

  # Same procedure of the WHS approch in order to compute the VaR
  sum = 0
  cont = -1
  while sum<=1-alpha:
    cont+=1
    sum+=loss_weights[cont]
  cont-=1

  VaR = loss[cont]

  return VaR*math.sqrt(riskMeasureTimeIntervalInYears*NumberOfDaysPerYears) #we already rescale the VaR up to 10 days

#Given function to compute a yearfrac using different time convention
def yearfrac(date1, date2, basis=0):
    """ Fraction of Year Between Dates.
    This function determines the fraction of a year occurring between two
    dates based on the number days between those dates using a specified
    day count basis. Based on MATLAB's yearfrac function.

    :param date1/date2:     values for dates (in pd.datetime format)
    :param basis:           2 for ACT/360
                            3 for ACT/365
                            6 for 30/360
    :return: year fraction
    """

    if basis == 2:
        return (date2-date1)/360
    elif basis == 3:
        return (date2-date1).days/365
    elif basis == 6:
        d2 = min(date2.day.any(), 30)
        d1 = min(date1.day.any(), 30)
        return (360*(date2.year-date1.year)+30*(date2.month-date1.month)+d2-d1) / 360
    else:
        print("Basis not recognised")
        return None

def cliquetOptionMonteCarlo(survP, bootstrap, stockPrice, L, volatility, notional, recovery, numIter):
  #defining parameters
  default_B = bootstrap.values[:,0]*survP
  nodef_B = bootstrap.values[:,0]
  fwd_B = np.array([ nodef_B[i]/nodef_B[i-1] for i in range(1,len(nodef_B))])
  yearFrac = yearfrac(bootstrap.index[0:-1],bootstrap.index[1:],3).to_numpy() #differenze tra le date
  deltaT = yearfrac(bootstrap.index[0],bootstrap.index[1:],3).to_numpy() #differenza da settlement date
  fwd_rate = - np.log(fwd_B) / yearFrac
  fwd_rate = np.insert(fwd_rate,0,bootstrap.values[0,1])


  mu = fwd_rate[0:7] #risk_free_rate

  MC = np.random.randn(7,numIter) #ogni colonna è simulazione montecarlo

  # computation of the stock dynamics with a GBM
  a = (mu-volatility**2/2)*yearFrac
  b = volatility*np.sqrt(yearFrac)
  lastPriceStock = stockPrice*np.ones((1,numIter))
  for i in range(7):
    new_price = lastPriceStock[-1,:]*np.exp( a[i] + b[i]*MC[i,:])
    lastPriceStock = np.vstack((lastPriceStock,new_price))


  # computation of the coupon using the simulated stock's dynamic
  coupon = np.array([ [max(L*lastPriceStock[i,j] - lastPriceStock[i-1,j],0) for j in range(numIter)] for i in range(1,8) ]) #cash flows
  defaultProb = np.array([survP[i-1] - survP[i] for i in range(1,len(survP))])

  discountedCoupon = np.array([coupon[:,i]*bootstrap.values[1:,0] for i in range(numIter)]).T #vedere per prodotto riga per colonna

  #explanation of the results in Assignment4-Group5_closedFormla.pdf
  discountedCoupon = np.flip(discountedCoupon,axis=0)
  sum = np.cumsum(discountedCoupon,axis=0)

  sum = np.flip(sum, axis=0)

  prices = [ np.dot(coupon[:,i],yearFrac*default_B[1:,]) for i in range(numIter)] + [recovery*np.dot(defaultProb,sum[:,i]) for i in range(numIter)]

  MCPrice = np.sum(prices)/numIter
  return MCPrice