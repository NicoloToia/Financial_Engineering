# !!! why the fuck is this a function????

def getTickers(df, names):
    """
    The function retrieves the tickers for the given names from the provided dataframe.

    INPUTS:
        df:     The dataframe containing the tickers and names
        names:  The names for which to retrieve the tickers

    OUTPUTS:
        tickers:  The list of tickers corresponding to the provided names
    """
    return df[df['Name'].isin(names)]['Ticker'].tolist()
