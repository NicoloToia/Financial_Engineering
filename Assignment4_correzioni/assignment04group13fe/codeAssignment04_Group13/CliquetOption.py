import numpy as np
from FE_Library import yearfrac
from scipy.stats import norm
import datetime


def CliquetOption(yearlySurvivalProb, dates, participationCoeff, notional, today, maturity, DiscountCurve, volatility, recovery):
    """
    The function calculates the price of a Cliquet option using the closed formula, which uses the Black's formula for a Call option

    INPUTS:
        yearlySurvivalProb:  The yearly survival probability of the counterparty
        dates:               The dates of the payoffs of the option
        participationCoeff:  The participation coefficient of the option
        notional:            The notional amount of the option
        today:               The current date
        maturity:            The maturity date of the option
        DiscountCurve:       The discount curve, with dates and corresponding discount factors
        volatility:          The volatility of the underlying asset
        recovery:            The recovery rate of the option

    OUTPUTS:
        price:  The price of the Cliquet option
    """
    yearlySurvivalProb = np.insert(yearlySurvivalProb, 0, 1)  # survival probability for today (P=1)
    dates = np.insert(dates, 0, today)
    discountFactors = np.interp(dates, DiscountCurve[:, 0], DiscountCurve[:, 1])  # obtain discount factors for dates from the discount curve
    dates_datetime = [datetime.datetime.fromordinal(int(date)) + datetime.timedelta(days=int(date % 1)) - datetime.timedelta(days=366) for date in dates]  # change the date format
    delta_t = [yearfrac(dates_datetime[i], dates_datetime[i + 1], 6) for i in range(maturity)]  # yearfracs between payment dates
    delta_t = np.array(delta_t)

    # compute via Black's formula the price of a call option with: initial value = L*B(0, t_n)/B(0, t_n-1), strike = 1, volatility
    d1 = ((np.log(participationCoeff * discountFactors[:maturity] / discountFactors[1:maturity+1]) + volatility**2 * delta_t / 2)
          / (volatility * np.sqrt(delta_t)))
    d2 = ((np.log(participationCoeff * discountFactors[:maturity] / discountFactors[1:maturity+1]) - volatility**2 * delta_t / 2)
          / (volatility * np.sqrt(delta_t)))
    c = (participationCoeff * norm.cdf(d1) - discountFactors[1:maturity+1] / discountFactors[:maturity] * norm.cdf(d2))

    price = np.sum(notional * c * ((delta_t * discountFactors[1:maturity+1] * yearlySurvivalProb[1:maturity+1]) +
             recovery * discountFactors[1:maturity+1] * (1 - yearlySurvivalProb[1:maturity+1])))

    return price
