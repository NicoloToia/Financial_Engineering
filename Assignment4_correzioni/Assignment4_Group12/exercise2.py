import pandas as pd
# import numpy as np
from assignment4functions import price_to_return
from assignment4functions import FullMonteCarloVaR
from assignment4functions import DeltaNormalVaR
from assignment4functions import SliceDataFromStartDate
from FE_Library import yearfrac


def runExercise2():
    """
    This function executes the exercise 2
    """
    # load EUROSTOXX_Dataset
    Dataset = pd.read_csv("EUROSTOXX50_Dataset.csv")

    # name of columns useful for Portfolio 2
    columns_name_2 = ["Date", "BMWG.DE"]

    # reduced dataset for Portfolio 2
    Dataset_2 = Dataset[columns_name_2]

    # Value date of VaR
    valueDate = '2017-01-16'
    date1 = pd.to_datetime(valueDate)

    # Maturity date of the call
    maturityDate = '2017-04-18'
    date2 = pd.to_datetime(maturityDate)

    # Slice the dataset from value date going backwards 2 years (for the WHS after)
    years = 2
    days_for_year = 365
    time_frame = years*days_for_year + 1
    Dataset_2 = SliceDataFromStartDate(Dataset_2, valueDate, time_frame)
    # as log returns I should take the ones of the last 2 years, going backwards from start date: 16.01.2015

    # compute returns for Portfolio 2: all daily returns until the value date
    logReturns = price_to_return(Dataset_2.copy())

    # Portfolio Value
    portfolioValue = 1186680

    # number of BMW stocks of the Portfolio 2
    stockPrice = Dataset_2.iloc[-1, 1]
    numberOfShares = portfolioValue/stockPrice

    # number of calls
    numberOfCalls = numberOfShares

    # Black formula for Call: parameters
    strike = 25
    timeToMaturityInYears = yearfrac(date1, date2, 3)
    rate = 0.005
    dividend = 0.031
    volatility = 0.154

    # Full MC VaR parameters
    NumberOfDaysPerYears = 256
    riskMeasureTimeIntervalInYears = 10/NumberOfDaysPerYears
    alpha = 0.95
    lambdaWHS = 0.95

    print('VaR with Full MC:')
    FullMCVaR = FullMonteCarloVaR(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend,
                                  volatility, timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha,
                                  NumberOfDaysPerYears, lambdaWHS)

    print('VaR with Delta-normal method:')
    DeltaNormal_VaR = DeltaNormalVaR(logReturns, numberOfShares, numberOfCalls, stockPrice, strike, rate, dividend,
                                     volatility, timeToMaturityInYears, riskMeasureTimeIntervalInYears, alpha,
                                     NumberOfDaysPerYears, lambdaWHS)

    return