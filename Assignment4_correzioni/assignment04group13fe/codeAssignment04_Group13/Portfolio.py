import numpy as np
from scipy.stats import norm
from scipy.stats import t

# Importing custom classes for different risk assessment and analysis methods
from HSMeasurements import HSMeasurements
from bootstrapStatistical import bootstrapStatistical
from WHSMeasurements import WHSMeasurements
from FullMontecarloVaR import FullMontecarloVaR
from PrincCompAnalysis import PrincCompAnalysis
from plausibilityCheck import plausibilityCheck
from DeltaNormalVaR import DeltaNormalVaR

# Class for European options modeling
from EuropeanOption import EuropeanOption


class Portfolio:
    """
        Represents a portfolio of assets, including both traditional assets and options. It allows for the
        computation of various risk measures (e.g., VaR, ES) using different methods. The portfolio can be constructed
        using either asset weights or the number of shares for each asset. Options can be added to the portfolio,
        and their impact on risk measures can be assessed only through some specific methods

        Attributes:
            assets (list):              List of asset symbols in the portfolio.
            dfAssetPrices (DataFrame):  Immutable historical price data for all assets, indexed by date. Used as a
                                        reference for analysis periods, adjustable via 'setDates'.
            options (dict):             Dictionary holding details of options added to the portfolio.
            numberOfDaysPerYears (int): Constant defining the number of trading days in a year used for calculations.
            portfolioValue (float):     The total current market value of the portfolio
            weights (numpy.array):      The proportion of each asset in the portfolio, calculated from either direct
                                        input or derived from the current market values, used in assessing the asset's
                                        contribution to the portfolio's overall risk and return.
            nShares (numpy.array):      The number of shares of each asset in the portfolio, set directly or calculated
                                        to reflect the portfolio's composition in terms of actual shares held.
            startDate (datetime):       The start date for the analysis period of the portfolio, set through 'setDates',
                                        marking the beginning of the historical data considered for calculations.
            endDate (datetime):         The end date for the analysis period, set through 'setDates', indicating the
                                        end of the period over which the portfolio is analyzed.
            prices (DataFrame):         A subset of 'dfAssetPrices', filtered to include only the data between
                                        'startDate' and 'endDate', used in all subsequent calculations requiring asset
                                        price data.
            returns (DataFrame):        Daily percentage returns of the assets, derived from 'prices'.
            logReturns (DataFrame):     Daily logarithmic returns of the assets, derived from 'prices'.
            meanReturns (Series):       Average daily percentage returns of each asset over the analysis period.
            covReturns (DataFrame):     Covariance matrix of the daily returns among all assets of interest.

        Methods:
            setPortfolioValue(self, portfolioValue):
                    Sets the total market value of the portfolio.
            updateNShares(self):
                    Updates the number of shares for each asset based on the current portfolio value and market prices.
            addEuOption(self, nContracts, optionType, underlying, strikePrice, riskFreeRate, dividend, maturity, volatility):
                    Adds a European option to the portfolio, specifying its details.

            VaR(self, alpha=0.99, riskMeasureTimeIntervalInDay=1, type='t', dof=1):
                    Calculates the Value at Risk (VaR) for the portfolio using either a normal distribution or a t-distribution.
            ES(self, alpha=0.99, riskMeasureTimeIntervalInDay=1, type='t', nu=4):
                    Computes the Expected Shortfall (ES), providing a more comprehensive risk measure than VaR.
            HSM(self, alpha=0.99, riskMeasureTimeIntervalInDay=1):
                    Utilizes historical simulation to calculate VaR and ES.
            bootstrapStat(self, numberOfSamplesToBootstrap=100, alpha=0.99, riskMeasureTimeIntervalInDay=1):
                    Applies the bootstrap method to estimate VaR and ES.
            WHSM(self, alpha=0.99, lambda_=1, riskMeasureTimeIntervalInDay=1):
                    Implements weighted historical simulation for risk measurement.
            PCA(self, alpha=0.99, numberOfPrincipalComponents=1, H=1):
                    Uses Principal Component Analysis (PCA) to assess risk, considering a specified number of principal components.
            pCheck(self, alpha=0.99, riskMeasureTimeIntervalInDay=1):
                    Conducts a plausibility check on the portfolio's returns.
            FullMC(self, UeOption, riskMeasureTimeIntervalInYears, alpha=0.95, NumberOfDaysPerYears=256):
                    Executes a full Monte Carlo simulation to assess the risk of the portfolio, including options.
            DNVaR(self, UeOption, riskMeasureTimeIntervalInYears, alpha=0.95, NumberOfDaysPerYears=256):
                    Determines the Delta-Normal Value at Risk, incorporating options into the risk assessment.
            printPortfolioDetails(self, method, VaR, ES=None):
                    Prints the details of the portfolio, including risk measures and composition.

        Private methods:
            __init__(self, assets, dfPrices, weightsORnShares, portfolioValue=None):
                    Initializes the portfolio with a list of assets, their prices, and either their weights or number of
                     shares. Optionally, the total portfolio value can be specified.
            __weightConstructor(self, weights, portfolioValue):
                    Private method for constructing the portfolio based on specified weights and total portfolio value.
            __nSharesConstructor(self, nShares):
                    Private method for constructing the portfolio based on the number of shares for each asset.
            __computeReturns(self):
                    Computes the daily percentage returns of the portfolio's assets.
            __computeLogReturns(self):
                    Calculates the log returns for more precise risk and performance analysis.
            __addOptions(self, ID, nContracts, maturity, derivative):
                    Private method to incorporate an option into the portfolio based on given parameters.setDates(self, start, end=0):
                    Sets the analysis period for the portfolio using start and end dates.
        """
    numberOfDaysPerYears = 256

    def __init__(self, assets, dfPrices, weightsORnShares, portfolioValue = None):
        """
        Initializes the Portfolio object with a list of assets, a DataFrame of asset prices, and either weights or
        the number of shares for each asset.
        Note: the total portfolio value can be specified for weight-based construction.

        IMPORTANT:
        'dfAssetPrices' is set with historical price data for these assets and is intended to remain immutable for
        direct modifications; adjustments to the portfolio's analysis period are made through 'setDates'.

        Parameters:
            assets (list): List of asset symbols.
            dfPrices (DataFrame): Immutable historical price data for all assets, intended not for modification but as
                                  a reference for specifying analysis periods through 'setDates'.
            weightsORnShares (list/array): Asset weights or number of shares, depending on the construction method.
            portfolioValue (float, optional): Total value of the portfolio, required if using weights to construct the portfolio.
        """
        self.assets = assets # Symbols for the portfolio's assets
        # Price data for the specified assets
        # This attribute is NOT intended to be modified: it is a shallow copy.
        self.dfAssetPrices = dfPrices[assets]
        self.options = {} # Initializes an empty dictionary for options

        # Sets the default start and end dates for price data analysis based on the DataFrame's date index
        # The method also inizialize the following variables: startDate, endDate, prices, meanReturns, covReturns
        self.setDates(self.dfAssetPrices.index[0], self.dfAssetPrices.index[-1])

        # Determines the construction path based on whether a total portfolio value is provided
        if portfolioValue is not None:
            self.__weightConstructor(weightsORnShares, portfolioValue)
        else:
            self.__nSharesConstructor(nShares=weightsORnShares)

    def __weightConstructor(self, weights, portfolioValue):
        """
        A private method to initialize the portfolio based on asset weights and total portfolio value. It calculates
        the number of shares for each asset based on these weights and the latest asset prices.

        Parameters:
            weights (array/list): The weights of each asset in the portfolio.
            portfolioValue (float): The total value of the portfolio.
        """
        self.setPortfolioValue(portfolioValue)
        self.weights = weights

        # Calculates the number of shares for each asset based on weights and portfolio value
        self.nShares = weights*portfolioValue/self.prices.loc[self.endDate]

    def __nSharesConstructor(self, nShares):
        """
        A private method to initialize the portfolio using the actual number of shares for each asset. It calculates
        the total value and the weights of the portfolio based on these shares and the latest asset prices.

        Parameters:
            nShares (array/list): The number of shares for each asset in the portfolio.
        """
        self.nShares = nShares

        # Calculates the total portfolio value based on the number of shares and latest asset prices
        self.setPortfolioValue( sum(nShares*(self.prices.loc[self.endDate])))

        # Calculates asset weights based on the number of shares
        self.weights = np.array(nShares)/sum(nShares)

    def __computeReturns(self):
        """
        A private method to calculate daily returns from price data. It updates the portfolio object with a DataFrame
        of these returns.
        """
        self.returns = self.prices.pct_change() # Computes percentage change for daily returns
        self.returns.dropna(inplace=True)       # Removes any NaN values from the DataFrame

    def __computeLogReturns(self):
        """
         A private method to calculate daily logarithmic returns from price data. It updates the portfolio object
         with a DataFrame of these log returns.
         """
        self.returns = np.log(self.prices / self.prices.shift(1)) # Computes log returns
        # !!!: this could be wrong
        self.returns.dropna(inplace=True)    # Removes any NaN values from the DataFrame

    def setDates(self, start, end=0):
        """
        Sets the start and end dates for the period over which the portfolio is analyzed. It then filters the
        portfolio's price data to this specific timeframe, ensuring that subsequent calculations only consider
        this period. Additionally, this method triggers the computation of daily and logarithmic returns based
        on the adjusted price data.

        Parameters:
            start (date-like): The start date for the analysis period. It should be compatible with pandas' date handling.
            end (date-like, optional): The end date for the analysis period. If not provided (or set to 0), the method
                                       uses the last available date in the price data as the end date.

        Note:
            - The method assumes that 'self.dfAssetPrices' contains a pandas DataFrame with a date index.
            - The method assumes that 'self.dfAssetPrices' does not contain any NaN values.
            - This method directly affects the 'self.prices' attribute, which is subsequently used for all
              financial calculations within the class.
        """
        # !!!: bruh why, just default to the self value
        # !!!: also relies on the ordering of dates
        if end == 0:  # Default end date handling: if not provided, use the last date from the asset prices DataFrame
            end = self.dfAssetPrices.index[-1]

        # Setting the internal start and end date attributes for future reference
        # !!!: qui hanno proprio cagato fuori dal vaso, DO NOT add or delete attributes dynamically
        self.startDate = start
        self.endDate = end

        # Filtering the asset prices DataFrame for the specified date range
        datesPrice = ((self.dfAssetPrices.index >= start) & (self.dfAssetPrices.index <= end))
        # !!!: here we are discarding information, what if we want to change dates again to a
        # !!!: wider range? We should just filter the prices when we need them or just pass already filtered prices
        self.prices = self.dfAssetPrices[datesPrice]

        # !!!: why do we recompute when dates are changed?
        # Initiating the computation of daily and logarithmic returns based on the filtered price data
        self.__computeReturns()     # Computes simple daily returns
        self.__computeLogReturns()  # Computes logarithmic returns

        # !!! why already? Should be done when needed (this is asking for a performance hit and trouble)
        # !!!: these should be computed only when needed and then cached
        self.meanReturns = self.returns.mean()  # Computes mean daily returns of the assets of interest
        self.covReturns = self.returns.cov()    # Computes covariance matrix of the returns of the assets of interest

        # !!!: I get this is called in the constructor but should really be a private method not
        # !!!: callable from outside

    def setPortfolioValue(self, portfolioValue):
        self.portfolioValue = portfolioValue

    def updateNShares(self):
        self.nShares = np.array((self.weights * self.portfolioValue / self.prices.loc[self.endDate]).values)

    def addEuOption(self, nContracts, optionType, underlying, strikePrice, riskFreeRate, dividend, maturity, volatility):
        """
        Adds a European option to the portfolio. This method calculates the time to maturity based on the portfolio's
        current end date and initializes a EuropeanOption object with the provided parameters. The option is then
        registered within the portfolio's options dictionary.

        Parameters:
            nContracts (int): The number of contracts for the European option.
            optionType (str): The type of the option ('call' or 'put').
            underlying (str): The symbol of the underlying asset for the option.
            strikePrice (float): The strike price of the option.
            riskFreeRate (float): The risk-free interest rate, assumed to be constant until the option's maturity.
            dividend (float): The dividend yield of the underlying asset.
            maturity (datetime): The maturity date of the option.
            volatility (float): The volatility of the underlying asset.

        Note:
            - This method assumes that the portfolio already contains the underlying asset and its price history.
            - The ID for the added option is constructed in a specific format for easy identification.
        """
        # !!!: should just pass a EuropeanOption object
        # !!!: why is this method? Why not use a more general version of addOption?

        # Fetching the last available price of the underlying asset from the portfolio's price data
        underlyingPrice = self.prices.iloc[-1][underlying].values[0]

        # Calculating time to maturity as a fraction of a year, based on the portfolio's defined number of trading days per year
        timeToMaturity = (maturity - self.endDate).days / self.numberOfDaysPerYears

        # Creating a new EuropeanOption object with the specified characteristics
        derivative = EuropeanOption(optionType, underlyingPrice, strikePrice, riskFreeRate, dividend, timeToMaturity, volatility)

        # Constructing a unique ID for the option
        ID = f"EU_{optionType}_{underlying[0]}"

        # Adding the option to the portfolio
        self.__addOptions(ID, nContracts, maturity, derivative)


    def __addOptions(self,ID, nContracts, maturity, derivative):
        """
        A private method to add the specified derivative (option) to the portfolio's options dictionary. It stores
        the option along with its details for further analysis and valuation.

        Parameters:
            ID (str): The unique identifier for the option.
            nContracts (int): The number of contracts of the option being added.
            maturity (datetime): The maturity date of the option.
            derivative (EuropeanOption): The EuropeanOption object representing the derivative to be added.

        Note:
            - This method is called internally by the addEuOption method.
            - The options dictionary is keyed by option IDs, which allows for efficient retrieval and management of options.
        """

        # !!!: options is not private, why is this method private?
        # !!!: better idea would be to just have a list rather than a dictionary

        # Adding or updating the option in the portfolio's options dictionary with its details
        self.options[ID] = {
            'nContracts': nContracts,
            'derivative': derivative,
            'maturity': maturity
            }

    def VaR(self, alpha=0.99, riskMeasureTimeIntervalInDay=1, type='normal', dof=1):
        """
        Computes the Value at Risk (VaR) for the traditional assets in the portfolio over a specified time interval and
        confidence level, using either a normal distribution or a t-distribution.

        Note: Options are not considered in this calculation.

        Parameters:
            alpha (float): Confidence level for VaR computation.
            riskMeasureTimeIntervalInDay (int): The time interval over which VaR is measured, in days.
            type (str): Type of distribution to use ('t' for t-distribution, any other value assumes normal distribution).
            dof (int): Degrees of freedom for the t-distribution (ignored if normal distribution is used).

        Returns:
            float: The calculated Value at Risk for the linear portfolio (it does not consider derivative)
        """

        #!!! using returns instead of logReturns
        # Calculate portfolio mean return and standard deviation
        mu = - self.weights.dot(self.meanReturns)
        std = np.sqrt( (self.weights.T.dot(self.covReturns)).dot(self.weights))

        # Compute VaR standard based on the distribution type
        if type == 't':
            VaR_std = t.ppf(alpha, df=dof)
        else: # if normal
            VaR_std = norm.ppf(alpha)

        # !!!: why does this return the VaR for the linear portfolio? Shouldn't it return the VaR for the whole portfolio?
        # !!!: bruh, you shit the bed on this, shoulda use the logreturns
        # Return the VaR value for the linear portfolio (it does not consider derivative)
        return self.portfolioValue * (riskMeasureTimeIntervalInDay * mu + np.sqrt(riskMeasureTimeIntervalInDay) * std * VaR_std)

    def ES(self, alpha=0.99, riskMeasureTimeIntervalInDay=1, type='normal', nu=4):
        """
        Computes the Expected Shortfall (ES) for the traditional assets in the portfolio over a specified time interval and
        confidence level, using a t-distribution or normal distribution.

        Note: Options are not considered in this calculation.

        Parameters:
            alpha (float): Confidence level for ES computation.
            riskMeasureTimeIntervalInDay (int): The time interval over which ES is measured, in days.
            type (str): Type of distribution to use ('t' for t-distribution, 'normal' for normal distribution).
            nu (int): Degrees of freedom for the t-distribution (ignored if normal distribution is used).

        Returns:
            float: The calculated Expected Shortfall for the linear portfolio (it does not consider derivative)
        """
        #!!! using returns instead of logReturns
        mu = - self.weights.dot(self.meanReturns)
        std = np.sqrt( (self.weights.T.dot(self.covReturns)).dot(self.weights))

        # Compute ES standard based on the distribution type
        if type == 't':
            t_inv_alpha = t.ppf(alpha, df=nu)
            phi_nu = t.pdf(t_inv_alpha, df=nu)
            ES_std = ((nu + t_inv_alpha ** 2) * phi_nu) / ((nu - 1) * (1 - alpha))
        else: # If normal
            # !!!: should be else if type == 'normal'
            z_alpha = norm.ppf(alpha)
            phi_z = norm.pdf(z_alpha)
            ES_std = phi_z / (1 - alpha)

        # !!!: error handling is missing

        # Return the ES value for the linear portfolio (it does not consider derivative)
        # !!!: should be D * mu + sqrt(D) * std * ES_std
        # ***: correct when D = 1
        return riskMeasureTimeIntervalInDay * self.portfolioValue * (mu + std * ES_std)

    def HSM(self, alpha=0.99, riskMeasureTimeIntervalInDay=1, print=False):
        """
        Computes the Value at Risk (VaR) and Expected Shortfall (ES) using the Historical Simulation Method (HSM).
        This approach uses the actual historical returns of the portfolio to estimate these risk measures.

        Note: Options are not considered in this calculation.

        Parameters:
            alpha (float): Confidence level.
            riskMeasureTimeIntervalInDay (int): The time period over which the risk is measured.

        Returns:
            !!! Why does this not consider derivatives? Shouldn't it return the VaR and ES for the whole portfolio?
            tuple: VaR, ES values for the linear portfolio (it does not consider derivative)
        """
        # !!! using returns instead of logReturns
        VaR, ES = HSMeasurements(self.returns, alpha, self.weights, self.portfolioValue, riskMeasureTimeIntervalInDay)
        if print:
            self.printPortfolioDetails('Historical Simulation Method', VaR, ES)
        return VaR, ES

    def bootstrapStat(self, numberOfSamplesToBootstrap=100, alpha=0.99, riskMeasureTimeIntervalInDay=1, print=False):
        """
        Computes VaR and ES by generating bootstrap samples from the portfolio's historical returns. This method
        provides a way to estimate risk measures from a potentially enhanced distribution of returns.

        Parameters:
            numberOfSamplesToBootstrap (int): The number of bootstrap samples to generate.
            alpha (float): Confidence level.
            riskMeasureTimeIntervalInDay (int): The time period over which the risk is measured.

        Returns:
            tuple: VaR, ES values for the portfolio.
        """
        # !!! using returns instead of logReturns
        samples = bootstrapStatistical(numberOfSamplesToBootstrap, self.returns)
        VaR, ES = HSMeasurements(samples, alpha,  self.weights, self.portfolioValue, riskMeasureTimeIntervalInDay)
        if print:
            self.printPortfolioDetails('Statistical Bootstrap', VaR, ES)
        return VaR, ES

    def WHSM(self, alpha=0.99, lambda_=1, riskMeasureTimeIntervalInDay=1, print=False):
        """
        Calculates VaR and ES using the Weighted Historical Simulation Method, which assigns different weights to
        historical returns, often emphasizing more recent data.

        Parameters:
            alpha (float): Confidence level.
            lambda_ (float): The decay factor to weight recent returns more heavily.
            riskMeasureTimeIntervalInDay (int): The time period over which the risk is measured.

        Returns:
            tuple: VaR, ES values for the portfolio.
        """
        #!!! using returns instead of logReturns
        VaR, ES = WHSMeasurements(self.returns, alpha, lambda_, self.weights, self.portfolioValue, riskMeasureTimeIntervalInDay)
        if print:
            self.printPortfolioDetails('Weighted Historical Simulation', VaR, ES)
        return VaR, ES

    def PCA(self, alpha=0.99, numberOfPrincipalComponents=1, H=1, print=False):
        """
        Utilizes Principal Component Analysis (PCA) to reduce the dimensionality of historical return data by focusing
        on the main sources of variance, for calculating VaR and ES.

        Parameters:
            alpha (float): Confidence level.
            numberOfPrincipalComponents (int): The number of principal components to consider.
            H (int): The horizon over which the risk is assessed, in days.

        Returns:
            tuple: VaR, ES values for the portfolio.
        """
        #!!! using returns instead of logReturns
        daysInYear = 256
        yearlyCovariance = daysInYear * self.covReturns
        yearlyMeanReturns = daysInYear * self.meanReturns
        VaR, ES = PrincCompAnalysis(yearlyCovariance, yearlyMeanReturns, self.weights, H, alpha, numberOfPrincipalComponents, self.portfolioValue)
        if print:
            self.printPortfolioDetails('PCA', VaR, ES)
        return VaR, ES

    def pCheck(self, alpha=0.99, riskMeasureTimeIntervalInDay=1, print = False):
        """
        Performs a plausibility check on the calculated VaR to assess its reasonableness given the historical
        return data of the portfolio.

        Parameters:
            alpha (float): Confidence level.
            riskMeasureTimeIntervalInDay (int): The time period over which the risk is measured.

        Returns:
            float: The plausibility-checked VaR for the portfolio.
        """
        #!!! using returns instead of logReturns
        VaR = plausibilityCheck(self.returns, self.weights, alpha, self.portfolioValue, riskMeasureTimeIntervalInDay)
        if print:
            self.printPortfolioDetails('Plausibility Check', VaR)
        return VaR

    def FullMC(self, UeOption, riskMeasureTimeIntervalInYears, alpha=0.95, NumberOfDaysPerYears=256, print=False):
        """
        Calculates VaR and ES using a full Monte Carlo simulation. This method is capable of incorporating the effects
        of options within the portfolio, making it suitable for complex portfolios with derivative instruments.

        Note: this method takes in account also European Option

        Parameters:
            UeOption (str): Identifier for the European option within the portfolio to be considered in the simulation.
            riskMeasureTimeIntervalInYears (float): The time period over which the risk is measured, in years.
            alpha (float): Confidence level.
            NumberOfDaysPerYears (int): The number of trading days in a year, used for scaling the simulation.

        Returns:
            float: The calculated VaR and possibly ES for the portfolio, incorporating the effects of the specified option.
        """
        numbersOfContracts = self.options[UeOption]['nContracts']
        option = self.options[UeOption]['derivative']

        VaR = FullMontecarloVaR(self.logReturns, self.nShares, numbersOfContracts, option.underlyingPrice, option.strike, option.riskFreeRate, option.dividend, option.volatility,
                      option.timeToMaturity, riskMeasureTimeIntervalInYears, alpha, NumberOfDaysPerYears)
        if print:
            self.printPortfolioDetails('Full MonteCarlo', VaR)
        return VaR


    def DNVaR(self, UeOption, riskMeasureTimeIntervalInYears, alpha=0.95, NumberOfDaysPerYears=256, print=False):
        """
        Calculates VaR using the Delta Normal method. This approach assumes normality of returns and linearity in
        portfolio exposures, including derivatives through their delta. It's quicker but less accurate than full
        Monte Carlo simulations.

        Note: this method takes in account also European Option

        Parameters:
            UeOption (str): Identifier for the European option within the portfolio to be considered in the calculation.
            riskMeasureTimeIntervalInYears (float): The time period over which the risk is measured, in years.
            alpha (float): Confidence level.
            NumberOfDaysPerYears (int): The number of trading days in a year, used for scaling the calculation.

        Returns:
            float: The calculated VaR for the portfolio, incorporating the effects of the specified option through delta approximation.
        """
        numbersOfContracts = self.options[UeOption]['nContracts']
        option = self.options[UeOption]['derivative']

        VaR = DeltaNormalVaR(self.logReturns, self.nShares, numbersOfContracts, option.underlyingPrice, option.strike, option.riskFreeRate, option.dividend, option.volatility,
                      option.timeToMaturity, riskMeasureTimeIntervalInYears, alpha, NumberOfDaysPerYears)
        if print:
            self.printPortfolioDetails('Delta Normal', VaR)
        return VaR

    def printPortfolioDetails(self, method, VaR, ES=None):
        """
        Prints the details of the portfolio, including the composition of shares and derivatives, the portfolio value,
        and the calculated Value at Risk (VaR) and Expected Shortfall (ES) based on the specified method.

        Parameters:
            method: The method used for calculating VaR and ES.
            VaR:    The calculated Value at Risk.
            ES:     The calculated Expected Shortfall (optional).

        Returns:
            None. This method prints the portfolio details to the console.
        """
        max_length = max(len(method), 15)  # 15 is the minimum width
        print("-------------------------------------------------------------")
        print("Portfolio composed by the following shares:")
        for ticker, weight in zip(self.assets, self.weights):
            print(f"\t{ticker:<15}\t{weight:<15.4f}")
        if self.options:
            print("Portfolio composed by the following derivatives:")
            for option_id, option_details in self.options.items():
                print(f"\t{option_id}")
        else:
            print(f"Portfolio value: {self.portfolioValue:>{max_length + 4}.4f} €")
        print(f"\n{method:<{max_length}} VaR: {VaR:>{15}.4f} €")
        if ES is not None:
            print(f"{method:<{max_length}} ES:  {ES:>{15}.4f} €")
