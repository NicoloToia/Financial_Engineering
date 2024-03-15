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
        d2 = min(date2.day, 30)
        d1 = min(date1.day, 30)
        return (360*(date2.year-date1.year)+30*(date2.month-date1.month)+d2-d1) / 360
    else:
        print("Basis not recognised")
        return None
