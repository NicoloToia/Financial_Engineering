# Review of theory
MATLAB is the core of the lessons. We will introduce a bit of new theory using State matrices and Markov chains.

## Notations
We use the notation of the book by Schonbrucher.

### Pure interest rate risk

$$ B(t,T) = \mathbb{E} \left[ e^{-\int_t^T r(s) \; ds} | \mathcal{F}_t \right] $$

$B$ is stochastic in the first variable but not in the second. It is a continuous family of stochastic functions, indexed by t.

#### Instantaneous forward rate

$$ f(t,T) = - \frac{\partial}{\partial T} \log B(t,T) $$

For us coupon Bond and forward rate will be synonymous.
Indeed we can write the following:

$$ B(t,T) = \exp \left( - \int_t^T f(t,u) \; du \right) $$

**Note**: here we are not selecting a model. We are simply characterising the forward rate as a function of time and maturity and the relation with discount factors.
We want the rates and discount factors, which are incredibly liquid, to be continuous.

**Attention**: the discount factor is based on a stochastic integral while the forward rate is a simple Rieman integral. In the discount factor we are explicitly writing out the expected value, while in the forward rate we are deriving a function of time $t$, since the expected value is a real valued function of time.

**Recall**: if $t=T$ the forward rate and instantaneous rate collapse to the same value. $$ f(t,t) = r(t) $$

## Market of Corporate bonds

In this market, we assume that a corporation can be unwilling or unable to pay back the bond. Here we employ stopping times and martingales to model the default of a corporation.

We introduce a random variable $\tau$, the stopping time, that models the default time. In particular we study the case of $t<\tau<T$. We will use the indicator function:

$$
\mathbb{I}(t) = \begin{cases}
    1 & \text{if } \tau > t \\
    0 & \text{if } \tau \leq t
\end{cases}
$$

The default is an absorbing state. Once the corporation defaults, it will not recover.

### Defaultable bond

In our case we model the defaultable zero coupon bond as follows:

$$
\bar{B}(t,T) = \mathbb{E} \left[
    \mathbb{I}(T) \cdot
    e^{-\int_t^T r(s) \; ds}
    | \mathcal{F}_t
\right]
$$

Like for the pure interest rate risk, we can define the instantaneous forward rate:

$$
\bar{f}(t,T) = - \frac{\partial}{\partial T} \log \bar{B}(t,T)
$$

and of course:

$$
\bar{B}(t,T) = \exp \left( - \int_t^T \bar{f}(t,u) \; du \right)
$$

## Relation between the two markets

By definition we have the following inequality:

$$
\bar{B}(t,T) \leq B(t,T)
\implies \bar{f}(t,T) \geq f(t,T)
$$

The second inequality is due to our assumption that the zero coupon bonds (either defaultable or not) be continuous. Thus we can split the integral into very small pieces and compare the two integrals term by term.

## Model

We are done with stylizing and laying the foundation of the theory. Now we pass to choosing and studying our model

#### Hypothesis
In particular we choose $r(t)$ and $\tau$ to be independent. This means we can factor the defaultable bond as follows:

$$
\bar{B}(t,T) = \mathbb{E} \left[
    \mathbb{I}(T) \cdot
    e^{-\int_t^T r(s) \; ds}
    | \mathcal{F}_t
\right]
= \mathbb{E} \left[
    \mathbb{I}(T)
    | \mathcal{F}_t
\right]
\cdot
\mathbb{E} \left[
    e^{-\int_t^T r(s) \; ds}
    | \mathcal{F}_t
\right]
$$

This is a simplification. It is a deliberate choice to model these two events as independent.

This does not reflect reality perfectly. But it highly simplifies many relations and leads to many closed form solutions.

Under this hypothesis we can define the following quantities:

$$
\begin{aligned}
    P(t,T) & = \mathbb{E} \left[ \mathbb{I}(T) | \mathcal{F}_t \right] \text{with } t \leq T \\
    \bar{B}(t,T) &= P(t,T) \cdot B(t,T)
\end{aligned}
$$

Under this model we can borrow some mathematical tools from the theory of Life Insurance.

### Instantaneous Hazard Rate

We define the following quantity:

$$
h(t,T) = \bar{f}(t,T) - f(t,T) \geq 0
$$

Furthermore:

$$
h(t,T) = - \frac{\partial}{\partial T} \log P(t,T)
$$

And this relation can be inverted to find $P(t,T)$:

$$
P(t,T) = \exp \left( - \int_t^T h(t,u) \; du \right)
$$

We can easily see that the difference between the two bonds is simply the hazard rate.

In other words, in a market or defaultable zero coupon bonds, in the limit case were the issuer walks away with no recovery, we can take a snapshot of the market (by Bootstrapping) to find the risk free term structure and afterwards and then writing everything in terms of the forward rate we can easily find the hazard rate.

Using the formalism of the forward rate opens the door to much more complex mathematical models. The hazard rate is the intensity of the risk of default.

Indeed, let us recall the definition of $\Lambda(t,T)$:

$$
P(t,T) = \exp \left( - \Lambda(t,T) \cdot (T-t) \right)
$$
If we choose $\Lambda(t,T) = \int_t^T h(t,u) \; du$ we can easily see that the two definitions are equivalent.

## Spot Zero Coupon Rate

We define the spot zero coupon rate as the quantity $y(t,T)$ such that:

$$
B(t,T) = \exp \left( - y(t,T) \cdot (T-t) \right)
$$

Thus:

$$
y(t,T) = - \frac{1}{T-t} \log B(t,T) = 
\frac{1}{T-t} \int_t^T f(t,u) \; du
$$

Again, we can think of the spot zero coupon rate as the intensity of the risk free zero coupon bond.

## Recovery Rate

In reality, defaults are not complete most of the time. We can recover part of our investment. To us, this will simply be a constant real value $\pi \in [0,1]$.

# Coupon Bearing Corporate Bond

Let us price a coupon bearing corporate bond. We will use the same model as before. The coupon has recovery rate $\pi$.

$$
\bar{C}(t) = 
    \sum_{n=1}^N \bar{C}_n \bar{B}(t, T_n) + \bar{B}(t, T_N)
    + \pi \cdot \sum_{n=1}^N \left[
        P(t, T_{n-1}) - P(t, T_n)
        \right] \cdot \bar{B}(t, T_n)
$$

with $P(t, T_0) = 1$ and $T_0 \leq t$. In other words, the we assume that the corporate will survive up to time $t$.

The first two terms are simply the present value of the coupons and the principal, using the defaultable bond. The third term is the present value of the eventual recovery (only of the principal) in case of default.

Unpaid coupons are not recovered.

**Note**: We are assuming default may only happen at coupon payment dates.

## Jarrow-Turnbull (1995) thumb rule
The default probability can be approximated byt the spread divided by 1 minus the recovery rate.

$$
P_{t,T} \approx \frac{s}{1-\pi}
$$

where $P_{t,T} = 1 - P(t,T)$ is the default probability and $s$ is the spread.

This way we get a family of default probability indexed by $\pi$.


### The Z-score

Usually, in practice the traders simply writes:

$$
\bar{C}(t) = \sum_{n=1}^N \bar{C}_n \hat{B}(t, T_n) + \hat{B}(t, T_N)
$$

where we define $\hat{B}(t,T)$ as follows:

$$
\hat{B}(t,T) = \exp \left( - \int_t^T \left[ f(t,u) + z(t,u) \right] \; du \right)
$$

where the quantity $z(t,T)$ is the zero spread and is chosen by the trader himself.

This is saying that there exists a spread such that the defaultable bond is equal to the default free bond plus the spread.

Naturally:

$$
0 \leq z(t,T) \leq h(t,T)
\implies
\bar{B}(t,T) \leq \hat{B}(t,T) \leq B(t,T)
$$

when $\pi = 0$ we have $\hat{B}(t,T) = B(t,T)$ and the two rates are the same.

**Note**: the z-spread is not a model. It is a deliberate choice of the trader. It is used to simplify the bootstrapping procedure.
We no longer have a family of defaultable bonds, but a single bond with a spread that is chosen by the trader.
$z(t,T)$ is a score, we rank bonds from most to least risky. $z(t,T)$ is a measure of the risk of the bond.

Here we have two sources of risk: the risk of default and the recovery rate.

For the market the default is a seamless event. We simply have an acceleraton of the coupon payments.

The z-score is a **SCORE** it tells us nothing about the underlying stochastic variables.

# Time-homogeneous Markov Chain
We suppose the states to be granular. In time the, we can move from one state to another. The probability is $q_{ij}$.
$i$ is the current state and $j$ is the next state.

We suppose that the probabilities stay costant in time.

# Rating Transition Matrix

This matrix gives us the probability of a bond to transition from one rating to another. It is a time-homogeneous Markov Chain.

We suppose the 

# Assignment

To reconstruct the hazard rate curve we must follow a common convention.
We have to make the curve piece-wise constant.

The z-spread is the parallel shift in the forward rate we need for the corporate bond to match the risk free bond.
 
If we have a curve of istantaneous forward rates, we can find the z-spread for each bond by shifting the curve by the hazard rate.
The z-spread is constant for each bond.

The z-spread is expressed in basis points. We must multiply by 10000 to get a number.

If we keep the hazard rate to be a constant, it can only have two values: one for the Investment Grade and one for the High Yield.

$$
h(t) = \begin{cases}
    h_{IG} & \text{if } R(t) = IG \\
    h_{HY} & \text{if } R(t) = HY
\end{cases}
$$

Now the hazard rate is a stochastic variable with only two variables. There derives a time-homogeneous Markov Chain.

We try to compute the transition matrix.

$$
\begin{bmatrix}
    q_{IG,IG} & q_{IG,HY} & q_{Default} \\
    q_{HY,IG} & q_{HY,HY} & q_{Default} \\
    0 & 0 & 1
\end{bmatrix}
$$