# Financial Engineering Projects - Group 16
This repository contains a series of assignments for the Financial Engineering course (2023-2024) developed by Group 16. The projects cover a wide range of topics related to financial derivatives, pricing models, and numerical methods. Below is a summary of the key topics covered in each assignment.

### Assignment 1: Option Pricing and Error Analysis

- European Call Option Pricing: Computed using Black’s formula, CRR binomial tree, and Monte Carlo (MC) methods.
- Error Rescaling: Evaluated the error in pricing using different methods and found optimal parameters for accuracy.
- Exotic Options: Pricing of a Knock-In Call option using closed-form solutions, CRR, and MC methods.
- Vega Sensitivity: Analyzed Vega for different options, comparing numerical and analytical results.
- Bermudan Option Pricing: Applied the CRR method to price a Bermudan option, and studied its behavior compared to vanilla European options.
- Monte Carlo Enhancements: Used antithetic variables to reduce variance in Monte Carlo simulations.
  
### Assignment 2: Yield Curves and Sensitivities

- Bootstrap Yield Curve Construction: Developed a discount factor curve and zero rates using market data, including deposits, futures, and swaps.
- Sensitivity Analysis: Computed DV01, BPV, and duration for interest rate swaps, analyzing the exposure of a portfolio to interest rate shifts.
- Theoretical Exercises: Derived bond pricing models and explored the Garman-Kohlhagen formula for European call options with time-dependent rates and volatility.
  
### Assignment 3: Asset Swaps, CDS, and Python Time-Series Analysis

- Asset Swap Spread: Calculated the Asset Swap Spread over Euribor 3m using market data.
- CDS Bootstrap: Constructed a CDS curve for various obligors (e.g., Intesa Sanpaolo) using bootstrapping and spline interpolation.
- First-to-Default Pricing: Priced a First-to-Default (FtD) contract between two obligors using the Li model and Monte Carlo simulations.
- Python Time-Series Analysis: Implemented basic time-series analysis in Python, including plotting log-returns and performing regression analysis.

### Assignment 4: Value at Risk (VaR) and Expected Shortfall (ES)
- Variance-Covariance Method: Computed VaR and ES at a 99% confidence level using a t-distribution for a portfolio of stocks.
- Historical Simulation and Bootstrap: Evaluated VaR and ES using historical simulation methods and explored their accuracy.
- Principal Component Analysis (PCA): Applied PCA to reduce dimensionality in VaR estimation.
- Monte Carlo VaR: Implemented full Monte Carlo simulation and Delta-Normal approaches for VaR analysis.
- Cliquet Option Pricing: Explored the pricing of a Cliquet option in the presence of counterparty risk​.

### Assignment 5: Advanced Derivative Pricing and Volatility Surface Calibration
- Certificate Pricing: Priced a certificate based on a basket of ENEL and AXA stocks using Monte Carlo simulations.
- Digital Option Pricing: Compared different methods (Black-Scholes, implied volatility) for pricing a digital option.
- Lewis Formula and FFT: Applied the Lewis formula and Fast Fourier Transform (FFT) for efficient pricing of call options.
- Volatility Surface Calibration: Calibrated a volatility surface using a mean-variance mixture model and analyzed the skew and vega of options​.

### Assignment 6: Interest Rate Risk and Hedging
- Bootstrap Market Discounts: Extended the bootstrapping technique to model discount factors and zero rates for interest rates beyond 12 years.
- Caplet Pricing and Spot Volatility Calibration: Calibrated spot volatilities from market cap prices using the Bachelier formula.
- Upfront Payment Calculation: Priced upfront payments for structured interest rate products.
- Delta and Vega Sensitivities: Computed delta and vega sensitivities for various buckets of maturities and strikes, and developed a hedging strategy using swaps and caps​.

### Assignment 7: Bermudan Swaptions and Certificate Pricing
- Bermudan Swaption Pricing: Priced a Bermudan swaption using the Hull-White model, exploring both closed-form solutions and tree-based methods.
- Certificate Pricing via NIG Model: Priced structured bonds using the Normal Inverse Gaussian (NIG) model, Fast Fourier Transform (FFT), and other methods (Variance Gamma and Monte Carlo).
- Black Model Adjustments: Evaluated the impact of digital risk adjustments on the Black model for certificate pricing.

## Energy Price and Load Forecasting

### EPLF Assignment 1: Regularization Techniques in Energy Price Forecasting
- Lasso, Ridge, and Elastic Net Regression: Explored regularization techniques to improve forecasting accuracy for electricity prices.
- Feature Selection: Demonstrated the importance of feature selection in time-series forecasting using LASSO.
- Seasonality Analysis: Investigated the effects of seasonality on energy price predictions using the Elastic Net model.

### EPLF Assignment 2: Deep Neural Networks (DNN) Hyperparameter Tuning
- Hyperparameter Optimization: Employed Optuna to perform random search for optimizing learning rate and hidden layer size in a DNN.
- Loss Function Minimization: Compared different configurations of learning rate and hidden size, evaluating their impact on the loss function.
- DNN Performance: Analyzed the effect of hyperparameters on overfitting and generalization performance in a DNN​.

### EPLF Assignment 3: Distributional Neural Networks and Quantile Regression
- Quantile Regression Neural Networks: Applied quantile regression to capture boundaries of data distribution using pinball loss functions.
- Distributional Neural Networks (DNN): Implemented DNNs using Normal and Johnson's SU distributions to capture probabilistic forecasting, tuning models via pinball and log-likelihood functions.
- Model Comparison: Compared model performance using Pinball and Winkler scores, finding JSU-DNN to perform best overall​.

## Risk Managment

### RM Assignment 1: Hazard Rate and Z-Spread Calculation
- Hazard Rate Curve Bootstrapping: Constructed hazard rate curves for investment grade (IG) and high yield (HY) bonds using bootstrapping and solving for hazard rates using MATLAB.
- Z-Spread Calculation: Calculated the Z-spread for IG and HY bonds by applying a parallel shift to the zero-rate curve to equalize defaultable and risk-free bond prices.
- Market-Implied Transition Matrix: Developed a time-homogeneous transition matrix for rating migrations, using the hazard rate curves​.

### RM Assignment 2: Present Value and Credit VaR
- Present Value Calculation: Evaluated present value using Monte Carlo simulations, considering survival probabilities and default scenarios for Investment Grade (IG) and High Yield (HY) bonds.
- Credit VaR Simulation: Simulated Credit VaR under different correlation scenarios, examining the impact of default and migration risks.
- Concentration Risk Analysis: Assessed the impact of concentration risk by varying the number of obligors, demonstrating the importance of diversification.
