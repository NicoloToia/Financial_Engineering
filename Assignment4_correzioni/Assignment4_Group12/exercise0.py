import pandas as pd
import numpy as np
from assignment4functions import AnalyticalNormalMeasures
from assignment4functions import price_to_return
from assignment4functions import SliceDataFromStartDate


def runExercise0():
    """
    This function executes the exercise 0
    """
    # load EUROSTOXX_Dataset
    Dataset = pd.read_csv("EUROSTOXX50_Dataset.csv")

    # name of columns useful for exercise 0
    columns_name = ["Date", "ADSGn.DE", "ALVG.DE", "MUVGn.DE", "OREP.PA"]

    # reduced dataset with companies to analyze
    Dataset_0 = Dataset[columns_name]

    # reduced Dataset_0
    end_date = '2020-02-20'
    years = 5
    days_for_year = 365
    time_frame = years*days_for_year + 1  # +1 because there's a leap year
    Dataset_0 = SliceDataFromStartDate(Dataset_0.copy(), end_date, time_frame)

    # create "matrix" of log-returns
    log_returns = price_to_return(Dataset_0.copy())

    # portfolio weights
    weights = np.full(4, 0.25)

    # significance value
    alpha = 0.99

    # value of portfolio (€)
    portfolioValue = 1.5e7  # 15000000

    # estimation interval
    riskMeasureTimeIntervalInDay = 1

    ES_portfolio, VaR_portfolio = AnalyticalNormalMeasures(alpha, weights, portfolioValue, riskMeasureTimeIntervalInDay,
                                                           log_returns)

    return