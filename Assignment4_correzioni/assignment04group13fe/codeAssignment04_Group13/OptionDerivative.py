class OptionDerivative:
    def __init__(self, underlyingPrice, strike, riskFreeRate, dividend, timeToMaturity, volatility):
        """
        Initializes an OptionDerivative object with parameters common to various option types.

        Parameters:
            underlyingPrice (float): The current price of the underlying asset.
            strike (float): The strike price of the option.
            riskFreeRate (float): The risk-free interest rate, expressed as a decimal.
            dividend (float): The dividend yield of the underlying asset, expressed as a decimal.
            timeToMaturity (float): The time to maturity of the option, in years.
            volatility (float): The volatility of the underlying asset's returns, expressed as a decimal.
        """
        self.underlyingPrice = underlyingPrice  # Current price of the underlying asset
        self.strike = strike                    # Strike price of the option
        self.riskFreeRate = riskFreeRate        # Risk-free interest rate
        self.dividend = dividend                # Dividend yield of the underlying asset
        self.timeToMaturity = timeToMaturity    # Time to maturity of the option, in years
        self.volatility = volatility            # Volatility of the underlying asset's returns

    # !!!: actually proper OOP requires a pure virtual function
    # ***: should add decorator @abstractmethod
    def price(self):
        pass
    def delta(self):
        pass
