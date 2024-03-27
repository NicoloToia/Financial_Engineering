import pandas as pd
import numpy as np
from FE_Library import yearfrac
from assignment4functions import blsCall
import random
import scipy.stats as stats
import math


def runExercise3():
    """
    This function executes the exercise 3
    """

    # input data
    partCoeff = 0.99  # participation coefficient
    notional = 3e7  # notional
    settlementDate = '2008-02-19'  # today's date
    number_years = 7  # expiration years for option
    sigma = 0.20  # volatility
    S0 = 1  # initial stock price
    recovery = 0.4  # recovery value

    # load boostrap DF values
    nameFile1 = 'Bootstrap.csv'
    bootstrapData = pd.read_csv(nameFile1, sep=';')
    bootstrapData['Date'] = pd.to_datetime(bootstrapData['Date'])
    bootstrapData.rename(columns={'Discount Factor': 'DF'}, inplace=True)

    # load boostrap CDS values
    nameFile2 = 'CDS_bootstrap.csv'
    survProbsData = pd.read_csv(nameFile2, sep=';')
    survProbsData['Date'] = pd.to_datetime(survProbsData['Date'])

    # Create vector of payment dates
    paymentDates = np.array(['2008-02-19', '2009-02-19', '2010-02-19', '2011-02-21', '2012-02-20', '2013-02-19',
                             '2014-02-19', '2015-02-19'])
    paymentDates = pd.to_datetime(paymentDates)

    # Do we have the discount value for each of these dates?
    # If not, we should perform a linear interpolation on zero rates
    # Let's obtain the zero rates
    bootStartDate = bootstrapData.iloc[0, 0]  # value date of bootstrap
    bootSize = len(bootstrapData['Date']) - 1  # number of bootstrapped dates (excluding the first one)
    yearFracBoot = np.zeros(bootSize)  # initialize yearfrac vec

    # Calculating needed year fractions
    for i in range(bootSize):
        yearFracBoot[i] = yearfrac(bootStartDate, bootstrapData.iloc[i + 1, 0], 3)

    # Zero Rates definition
    zeroRates = - np.log(bootstrapData.iloc[1:, 1]) / yearFracBoot

    # Year fractions for payment dates
    numPayments = int(len(paymentDates) - 1)
    yearFracPayment = np.zeros(numPayments)
    for i in range(numPayments):
        yearFracPayment[i] = yearfrac(bootStartDate, paymentDates[i + 1], 3)

    # Interpolation to obtain the zero rates for all the payments dates
    zeroRatesPayment = np.interp(yearFracPayment, yearFracBoot, zeroRates)

    # Find the discount factors for payment dates
    discFactorsPayment = np.exp(-zeroRatesPayment * yearFracPayment)
    # we checked the results with the previous lab on matlab

    # Compute default probabilities (year by year)
    survProbs = survProbsData.iloc[:, 1].to_numpy()  # we will need it as numpy array afterward
    defaultProbs = np.add(np.array([1] + survProbs[:-1].tolist()), - survProbs)
    # yearly default probabilities(year 0 to 6)

    # Find the forward discount rates
    forwardDiscFact = discFactorsPayment / np.array([1] + discFactorsPayment[:-1].tolist())  # fwd discount factors
    # explanation of code:  discFactorsPayment[:-1].tolist() converts from np array to list,
    # so that we can append it to the element 1, using the +, creating a new list,
    # then turn it back to numpy; 1 is the DF B(0,0)

    deltaPayments = np.zeros(numPayments)
    for i in range(numPayments):
        deltaPayments[i] = yearfrac(paymentDates[i], paymentDates[i + 1], 3)  # yearfracs between payment dates

    forwardDiscRates = -np.log(forwardDiscFact) / deltaPayments  # forward discount rates

    #####################
    # WE HAVE OBTAINED ALL THE DATA WE WILL NEED TO USE AFTERWARDS
    #####################

    # Analytical price for the Cliquet option
    # We consider the Cliquet option as a sum of individual Forward Start Call options, each one lasting one year
    # and starting when the previous one expires

    # Price each individual Forward Start Call: for derivation and explanation of formula please consult our report
    # We already multiply by the notional
    callPricesVec = blsCall(S0, S0 / partCoeff, forwardDiscRates, 0, sigma, deltaPayments) * partCoeff * notional

    # Sum the single option prices, multiply by the notional, and we obtain the analytical price with NO default
    priceNoDefaultAN = np.sum(callPricesVec)
    print('Analytical price NON-defaultable: ', priceNoDefaultAN)

    # Create vector of yearly future CF, as before
    futureCashFlows = np.flip(np.cumsum(callPricesVec))  # future CF for year 0 to year 6

    # Price in defaultable case
    priceDefaultAN = np.sum(callPricesVec * survProbs) + recovery * np.dot(futureCashFlows, defaultProbs)
    print('Analytical price defaultable: ', priceDefaultAN)

    ####################
    # MC Price for the Cliquet option - 1st method: MC for time to default and then analytical formula

    # MC parameters
    nSim = int(1e4)
    random.seed(40)

    # draw random probabilities
    rand_U = np.random.uniform(0, 1, (nSim, 1))

    # Time to defaults in a smart way
    countNoDef = 0  # the price for these "paths" is always the same because we receive all coupons,
    # so we just need to know how many they are
    tauDef = np.array([])  # just to initialize the array

    for i in range(nSim):
        indexYear = np.where(survProbs <= rand_U[i])  # we find where the tau falls
        # I'm finding the index which corresponds to the first prob. with value smaller than this
        if len(indexYear[0]) == 0:  # when tau_ISP>maturity
            countNoDef += 1
        else:
            tauDef = np.append(tauDef, indexYear[0][0])  # append the reset date when default is "computed"

    # Evaluate the option value considering the different time to defaults, obtaining one option price for each sim
    # We use the analytical formula as in the other part of this exercise

    # No-default case
    priceSimsMC = np.full(countNoDef, priceNoDefaultAN)  # we need a vector to compute the std of MC price after

    # Default case: we have to see when our simulation defaults and consider the necessary cash flows in response
    for i in range(len(tauDef)):
        if tauDef[i] == 0:
            partialSum = futureCashFlows[0] * recovery  # defaults between 0 and 1, so we don't get any "pure coupon"
        else:
            partialSum = np.sum(callPricesVec[:int(tauDef[i])]) + recovery * futureCashFlows[int(tauDef[i])]
            # add contributes of "pure coupons" before default and then recovery part for future expected CF
        priceSimsMC = np.append(priceSimsMC, partialSum)  # appending the price computed for this simulation

    # Compute final MC price via the mean of the summed prices of all simulations
    priceDefaultMC = np.mean(priceSimsMC)
    print('MC price Defaultable 1st Method: ', priceDefaultMC)

    # Compute sample standard deviation
    stdVec = np.sum((np.add(priceSimsMC, - priceDefaultMC)) ** 2)
    stdVec = np.sqrt(stdVec / (nSim - 1))

    # Compute MC confidence interval
    normQuantile = stats.norm.ppf(q=0.975, loc=0, scale=1)  # quantile z(alpha/2) of standard normal
    lowBoundMC = priceDefaultMC - normQuantile * stdVec / math.sqrt(nSim)
    upBoundMC = priceDefaultMC + normQuantile * stdVec / math.sqrt(nSim)
    confIntMC = [lowBoundMC, upBoundMC]
    print('MC price Defaultable confidence interval 1st Method: ', confIntMC)

    ###################################
    ###################################
    # FROM HERE ON IT'S JUST SOME 2 OTHER METHODS TRIED, BUT WE WON'T REPORT THEIR RESULTS
    ###################################
    ###################################

    # MC Price for the Cliquet option - 2nd method: MC for underlying evolution and then analytical formula for defaults

    # Random realizations of std Gaussian r.v. to simulate underlying's GBM dynamics
    gaussRV = np.random.normal(loc=0, scale=1, size=(nSim, numPayments))  # numPayments is 7

    # Initializing matrix to store simulated underlying value in time
    underlyingDynamics = np.ones((nSim, numPayments + 1))
    underlyingDynamics[:, 0] = S0 * underlyingDynamics[:, 0]

    # Evaluating the underlying simulated dynamics: Geometric Brownian Motion
    for i in range(numPayments):
        underlyingDynamics[:, i + 1] = underlyingDynamics[:, i] * np.exp(
            (forwardDiscRates[i] - sigma ** 2 / 2) * deltaPayments[i]
            + sigma * np.sqrt(deltaPayments[i]) * gaussRV[:, i])

    # Now that we have the stock dynamics in time, let's evaluate the option's payoff,
    # year by year, one payoff for each simulation
    payoffMC = np.maximum(0, np.add(partCoeff * underlyingDynamics[:, 1:], - underlyingDynamics[:, :-1]))

    # Compute discounted payoff, one for each year
    # (one for each fwd start call constituting our Cliquet), including notional
    payoffDiscMC = np.mean(payoffMC, axis=0) * discFactorsPayment * notional
    # compute the mean payoff of each column (so one for each year) and then discount them by the corresponding DF

    # Compute sample standard deviation
    stdVec = np.zeros(len(payoffDiscMC))
    for i in range(len(payoffDiscMC)):
        stdVec[i] = np.sum((np.add(payoffMC[:, i], - np.mean(payoffMC, axis=0)[i])) ** 2)
    sampleStdMC = np.sqrt(stdVec / (nSim - 1)) * discFactorsPayment * notional
    # sample standard deviation computed as array that will then be summed (one std for each call price)
    # we have to multiply by the discount factors and the notional to be "financially" coherent with the payoff value

    # Evaluate final MC price of Cliquet option in case of no counterparty risk
    priceNoDefaultMC = np.sum(payoffDiscMC)
    # print('MC price NON-defaultable 2nd Method: ', priceNoDefaultMC)

    # Evaluate MC confidence interval for the Cliquet option price in case of no counterparty risk
    normQuantile = stats.norm.ppf(q=0.975, loc=0, scale=1)  # quantile z(alpha/2) of standard normal
    lowBoundNoDefMC = np.sum(np.add(payoffDiscMC, - normQuantile * sampleStdMC / math.sqrt(nSim)))
    upBoundNoDefMC = np.sum(np.add(payoffDiscMC, + normQuantile * sampleStdMC / math.sqrt(nSim)))
    confIntNoDefMC = [lowBoundNoDefMC, upBoundNoDefMC]
    # print('MC price NON-defaultable confidence interval 2nd Method: ', confIntNoDefMC)

    # Now we can evaluate the price of the option in case of possible defaults
    # How do we consider the recovery part?
    # We already have the strip of future default-free discounted cashflows, we just need to weight them
    # by the default probability (for each interval) and then multiply them by recovery value R
    # Example: if we default between t2 and t3, we will apply the recovery to all the future cf from that point,
    # so we are excluding the already received cash flows from the recovery part

    # To do this in a smart way, we will create a vector of future cash flows for each year
    # (gets smaller the more we go on)
    futureCashFlowsMC = np.flip(np.cumsum(payoffDiscMC))  # future CF for year 0 to year 6

    # Evaluate MC price of Cliquet option in presence of counterparty risk
    priceDefaultMC = np.dot(survProbs, payoffDiscMC) + recovery * np.dot(defaultProbs, futureCashFlowsMC)
    # print('MC price defaultable 2nd Method: ', priceDefaultMC)

    # Evaluate MC confidence interval of Cliquet option in presence of counterparty risk
    lowStdPayoffDiscMC = np.add(payoffDiscMC, - normQuantile * sampleStdMC / math.sqrt(nSim))
    upStdPayoffDiscMC = np.add(payoffDiscMC, + normQuantile * sampleStdMC / math.sqrt(nSim))

    lowStdFutureCF = np.flip(np.cumsum(lowStdPayoffDiscMC))
    upStdFutureCF = np.flip(np.cumsum(upStdPayoffDiscMC))

    lowBoundDefMC = np.dot(survProbs, lowStdPayoffDiscMC) + recovery * np.dot(defaultProbs, lowStdFutureCF)
    upBoundDefMC = np.dot(survProbs, upStdPayoffDiscMC) + recovery * np.dot(defaultProbs, upStdFutureCF)

    confIntDefMC = [lowBoundDefMC, upBoundDefMC]
    # print('MC price defaultable confidence interval 2nd Method: ', confIntDefMC)

    ######################
    # MC Price for the Cliquet option - 3rd method: MC for time to default and MC for underlying evolution
    # MC parameters
    nSim = int(1e4)

    # draw random probabilities
    rand_U = np.random.uniform(0, 1, (nSim, 1))

    # Time to defaults in a smart way
    countNoDef = 0  # the price for these "paths" is always the same because we receive all coupons,
    # so we just need to know how many they are
    tauDef = np.array([])  # just to initialize the array

    for i in range(nSim):

        indexYear = np.where(survProbs <= rand_U[i])  # we find where the tau falls
        # I'm finding the index which corresponds to the first prob. with value smaller than this

        if len(indexYear[0]) == 0:  # when tau_ISP>maturity
            countNoDef += 1
        else:
            tauDef = np.append(tauDef, indexYear[0][0])

    # Now we evaluate the option value considering the different time to defaults,
    # obtaining one option price for each simulation
    # We use the analytical formula as in the other part of this exercise

    # No-default case
    sumPriceSimsMC = priceNoDefaultMC * countNoDef

    # Default case: we have to see when our simulation defaults and consider the necessary cash flows in response
    for i in range(len(tauDef)):

        if tauDef[i] == 0:
            partialSum = futureCashFlowsMC[
                             0] * recovery  # defaults between 0 and 1, so we don't receive any "pure coupon"
        else:
            partialSum = np.sum(payoffDiscMC[:int(tauDef[i])]) + recovery * futureCashFlowsMC[int(tauDef[i])]
            # add contributes of "pure coupons" before default and then recovery part for future expected CF
        sumPriceSimsMC += partialSum

    # Compute final MC price via the mean of the summed prices of all simulations
    priceDefaultMC = sumPriceSimsMC / nSim
    # print('MC price Defaultable 3rd Method: ', priceDefaultMC)

    return
