# Probabilistic Forecasting Lab

## Introduction

We need our probabilistic forescasting to, of course, really capture the probability distribution of the data:

$$
\mathbb{P}(y_{n+1} \in \mathcal{C}(x_{n+1})) \approx 1 - \alpha
$$

therefore, if we claim to be capturing 80% of the data, we should be capturing 80% of the data and not less or more. This is the main goal of this lab.

Of course, we should also capture how the data is distributed.
For example, the data could have wider spread in certain regions of the input space, and we should be able to capture that too.
This could be done, for example, by using another distribution for the variance in the y space.

What we are trying to predict is the distribution of the output, given the input.

$$
\mathbb{P}(Y | X)
$$

This can be achieved by creating a single NN that outputs the parameters of the distribution of all outputs.
We could also have a NN for each hour of the day and estimate the distribution of the output of that hour.

## Quantile regression NN

We have the following equations:

$$
\begin{align*}

\ell_1 &= g(x_i \; W_1 + b_1) \\
\ell_2 &= g(\ell_1 \; W_2 + b_2) W_3 + b_3 \\

& W_1 \in R^{n_x \times n_{u_1}}, \quad W_2 \in R^{n_{u_1} \times n_{u_2}}\\
& W_3 \in R^{n_{u_2} \times H \cdot n_p}, \quad n_{u_1}, n_{u_2} \in \mathbb{Z}^+\\
& b_1 \in R^{n_{u_1}}, \quad b_2 \in R^{n_{u_2}}, \quad b_3 \in R^{H \cdot n_p}\\

\end{align*}

$$

In our case we try to approximate the deciles of our data's distribution.
In order to do that, we need an objective function that captures the quantiles of the data.
We use the pinball loss function:

$$

\begin{align*}

n_p = \# \Gamma (\text{number of quantiles, e.g. 10 for the declies}) \\

\sum_i \sum_h \sum_{\gamma} (y_i^h - \hat{q}_{\gamma}^h(x_i))^+ \gamma + (y_i^h - \hat{q}_{\gamma}^h(x_i))^- (1 - \gamma)

\end{align*}
$$

Where $y_i^h$ is the true value of the output at hour $h$, $x_i$ is the input, $\hat{q}_{\gamma}^h(x_i)$ is the predicted quantile of the output at hour $h$ and $\gamma$ is the quantile we are trying to predict.

## Distributional NN

We have the following equations:

$$
\begin{align*}

\ell_1 &= g(x_i \; W_1 + b_1) \\
\ell_2 &= g(\ell_1 \; W_2 + b_2) W_3 + b_3 \\

& W_1 \in R^{n_x \times n_{u_1}}, \quad W_2 \in R^{n_{u_1} \times n_{u_2}}\\
& W_3 \in R^{n_{u_2} \times H \cdot n_p}, \quad n_{u_1}, n_{u_2} \in \mathbb{Z}^+\\
& b_1 \in R^{n_{u_1}}, \quad b_2 \in R^{n_{u_2}}, \quad b_3 \in R^{H \cdot n_p}\\

\end{align*}

$$

For example, for Johnson's SU distribution, we have the following parameters:

$$
\begin{cases}

    \lambda_i^h = \ell_2^{[h]} \\
    \sigma_i^h = \epsilon + \gamma \text{Softplus}(\ell_2^{[H+h]}) \\
    \tau_i^h = 1 + \gamma \text{Softplus}(\ell_2^{[2H+h]}) \\
    \zeta_i^h = \ell_2^{[3H+h]} \\
    \text{Softplus}(x) = \log(1 + \exp(x))
\end{cases}
$$

Here instead of computing the quantiles we have as output the parameters of the distribution of the data.

Our objective function will now be a log-likelihood ratio test (Negative Log Likelihood).

Of course, since NN are just simple non-linear function approximators, we can use any one of them to approximate the distribution of the data.

The key point in the porbabilistic approach is the output layer.
What our model outputs is not the mean of the distribution, but the parameters of the distribution.

All other techniques (pruning,  tuning, ecc. ecc.) we have seen so far apply to this kind of models too.

## Coding

### Quantile Crossing

It can happen that when we are estimating  quantiles quite close to each other, they might cross each
other.
This can be a problem.
A very simple method to avoid is this is to simply sort the quantiles, thus switching the order of
the quantiles.

### Objective Function

We remark that whatever we do to build our NN (or other model) is stricly independent of the
objective function we are using.

### Output Layer for testing

With the point prediction our ouput was a vector of 24 elements, one for each hour of the day.
In the quantile regression it is now the number of quantiles.
In the distributional approach it is the number of parameters of the distribution for each hour of the day
(e.g. 2x24) or again the number of quantiles.

The type of model we are using change the size of our output layer.

### Exper_config.json

In the exper_config.json file we can set the method to use and the number of quantiles to regress.

# HW:

Experiment if the Johnson's SU distribution is better than the quantile regression or the Normal distribution.
Change the output size to fit the Johnson's SU distribution (4 parameters) and compare the results.

Implement the matrix reshaping for the Johnson's SU distribution.

Perform a random search on the hyperparameters of the Johnson's SU distribution.
Test a recalibration to use May 2017 as a test set.

## Facultative HW

Implement the Winkler score.
Preprocess the price by using the *arcsinh* function and experiment on the QR.
We must apply the inverse transformation to the predictions before computing the Winkler score.
