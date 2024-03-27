import pandas as pd
import numpy as np
from assignment4functions import price_to_return
from assignment4functions import HSMeasurements
from assignment4functions import bootstrapStatistical
from assignment4functions import WHSMeasurements
from assignment4functions import plausibilityCheck
from assignment4functions import SliceDataFromStartDate
from assignment4functions import PrincCompAnalysis


def runExercise1():
    """
    This function executes the exercise 1
    """
    # load EUROSTOXX_Dataset
    Dataset = pd.read_csv("EUROSTOXX50_Dataset.csv")

    print('Portfolio A')

    # name of columns useful for Portfolio A
    columns_name_A = ["Date", "AXAF.PA", "SASY.PA", "TTEF.PA", "VOWG_p.DE"]

    # reduced dataset for Portfolio A, taking only the needed columns
    Dataset_A = Dataset[columns_name_A]

    # reduced Dataset_A, taking only the needed rows
    end_date = '2019-03-20'
    years = 5
    days_for_year = 365
    time_frame = years * days_for_year + 1
    Dataset_A = SliceDataFromStartDate(Dataset_A.copy(), end_date, time_frame)

    # compute returns for Portfolio A
    returns_A = price_to_return(Dataset_A.copy())

    # computation of weights and the value of the portfolio
    sharesNumber = np.array([20e3, 20e3, 25e3, 10e3])  # number of shares bought of each stock - alphabetical order
    todayPrices = Dataset_A[Dataset_A['Date'] == end_date]  # take the prices stockPrice(t) at value date ( = end_date)
    todayPrices = todayPrices.iloc[:, 1:].to_numpy().reshape(len(sharesNumber))
    # convert it to numpy array and reshape as sharesNumber (4,)
    portfolioValue_A = (sharesNumber * todayPrices).sum()  # initial value of Portfolio
    weights = sharesNumber * todayPrices / portfolioValue_A  # FROZEN PORTFOLIO assumption

    # significance value
    alpha = 0.95
    # estimation interval
    riskMeasureTimeIntervalInDay = 1

    # Historical simulation
    print('Var and ES with historical simulation')
    ES_A_HS, VaR_A_HS = HSMeasurements(returns_A, alpha, weights, portfolioValue_A, riskMeasureTimeIntervalInDay)

    # Statistical bootstrap
    numberOfSamplesToBootstrap = 200
    samples = bootstrapStatistical(numberOfSamplesToBootstrap, returns_A)  # randomly selecting indexes
    returns_A_bootstrap = returns_A.iloc[samples, :]  # selecting the corresponding returns

    print('Var and ES with statistical bootstrap')
    ES_A_bootstrap, VaR_A_bootstrap = HSMeasurements(returns_A_bootstrap, alpha, weights, portfolioValue_A,
                                                     riskMeasureTimeIntervalInDay)

    # check on the order of magnitude
    VaR_HS_ptf1_check = plausibilityCheck(returns_A, weights, alpha, portfolioValue_A, riskMeasureTimeIntervalInDay)

    print('')
    ##############################################################################################################
    print('Portfolio B')

    # name of columns useful for Portfolio B
    columns_name_B = ["Date", "ADSGn.DE", "AIR.PA", "BBVA.MC", "BMWG.DE", "DTEGn.DE"]

    # reduced dataset for Portfolio B, taking only the needed columns
    Dataset_B = Dataset[columns_name_B]

    # reduced Dataset_C, taking only the needed rows
    end_date = '2019-03-20'
    years = 5
    days_for_year = 365
    time_frame = years * days_for_year + 1
    Dataset_B = SliceDataFromStartDate(Dataset_B.copy(), end_date, time_frame)

    # compute returns for Portfolio B
    returns_B = price_to_return(Dataset_B.copy())

    # Portfolio B weights (since we are investing 20 cents in each stock)
    portfolioValue_B = 1
    weights_B = np.full(5, 0.20)  # FROZEN PORTFOLIO assumption

    # significance value
    alpha_B = 0.95
    # historical index
    lambda_B = 0.95
    # estimation interval
    riskMeasureTimeIntervalInDay_B = 1

    # evaluate VaR and ES for Portfolio B with WHS
    print('Var and ES with weighted historical simulation')
    ES_B_WHS, VaR_B_WHS = WHSMeasurements(returns_B, alpha_B, lambda_B, weights_B, portfolioValue_B,
                                          riskMeasureTimeIntervalInDay_B)

    # check on the order of magnitude
    VaR_WHS_ptf2_check = plausibilityCheck(returns_B, weights_B, alpha_B, portfolioValue_B,
                                           riskMeasureTimeIntervalInDay_B)

    print('')
    ###########################################################################################################
    print('Portfolio C')

    # load indexes.csv
    indexes = pd.read_csv('_indexes.csv')
    # name of columns useful for Portfolio C from indexes.csv
    columns_name_C = indexes['Ticker'].tolist()[:3] + indexes['Ticker'].tolist()[4:19]
    # add "Date"
    columns_name_C.insert(0, "Date")

    # reduced dataset for Portfolio C, taking only the needed columns
    Dataset_C = Dataset[columns_name_C]

    # reduced Dataset_C, taking only the needed rows
    end_date = '2019-03-20'
    years = 5
    days_for_year = 365
    time_frame = years * days_for_year + 1
    Dataset_C = SliceDataFromStartDate(Dataset_C.copy(), end_date, time_frame)

    # compute returns for Portfolio
    returns_C = price_to_return(Dataset_C.copy())

    # compute covariance matrix YEARLY
    cov_matrix = pd.DataFrame.cov(returns_C.iloc[:, 1:])  # DAILY
    buss_days_year = 256
    cov_matrix = cov_matrix * buss_days_year

    # compute mean YEARLY
    mean_vector = returns_C.iloc[:, 1:].mean(axis=0)  # done by columns -> one mean for each company
    mean_vector = mean_vector * buss_days_year

    # Portfolio C weights (since we are investing 1 euro in total)
    weights_C = np.full(18, 1 / 18)
    portfolioValue_C = 1

    # significance value
    alpha_C = 0.95

    # days YEARLY
    H = 10 / buss_days_year

    # initializing vectors that will contain the values of ES and VaR for different number of principal components
    ES_C_PCA = np.zeros(5)
    VaR_C_PCA = np.zeros(5)

    print('Var and ES with principal component analysis')
    # number of principal components goes from 1 to 5
    for number_pc in range(1, 6):
        ES_C_PCA[number_pc-1], VaR_C_PCA[number_pc-1] = PrincCompAnalysis(cov_matrix, mean_vector, weights_C, H,
                                                                          alpha_C, number_pc, portfolioValue_C)

    print('VaR with number of components from 1 to 5:', VaR_C_PCA)
    print('ES with number of components from 1 to 5:', ES_C_PCA)

    # check on the order of magnitude on the VaR
    VaR_WHS_ptf2_check = plausibilityCheck(returns_C, weights_C, alpha_C, portfolioValue_C, H*buss_days_year)

    return