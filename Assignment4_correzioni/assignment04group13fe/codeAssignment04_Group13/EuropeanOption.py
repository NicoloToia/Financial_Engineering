from OptionDerivative import OptionDerivative
from blsprice import blsprice
from blsprice import blsdelta

class EuropeanOption(OptionDerivative):
    def __init__(self, optionType, *args, **kwargs):
        """
        Initializes a EuropeanOption object as a specific type of OptionDerivative, incorporating all the
        parameters from the parent class and adding the option type (call or put).

        Parameters:
            optionType (str): Specifies the type of the option - 'call' for call options or 'put' for put options.
        """
        super().__init__(*args, **kwargs)
        self.optionType = optionType  # Specifies the option type ('call' or 'put')

    def price(self):
        """
        Calculates the price of the European option based on the Black-Scholes formula.

        Returns:
            float: The calculated price of the European option.
        """

        if self.optionType == 'call':
            flag = 1    # Indicates a call option to the pricing function
        else:  # put
            flag = -1   # Indicates a put option to the pricing function

        # Call the blsprice function with parameters inherited from OptionDerivative and the option type flag
        return blsprice(self.underlyingPrice, self.strike, self.riskFreeRate, self.dividend, self.volatility,
                        self.timeToMaturity, flag)

    # !!!: missing error handling for invalid optionType, should raise ValueError

    def delta(self):
        """
        Calculates the delta of the European option
        This implementation assumes a similar approach for both call and put options, based on the Black-Scholes model.

        Returns:
            float: The calculated delta of the European option.
        """

        return blsdelta(self.underlyingPrice, self.strike, self.dividend, self.volatility, self.timeToMaturity)
    
    # !!!: same as above
