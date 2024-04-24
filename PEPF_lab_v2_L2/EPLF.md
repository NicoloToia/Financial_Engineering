# Intro to DNN

## Gradient Descent

Our loss function is:

$$
R(\theta) = \frac{1}{2} \sum_{i=1}^{n} (y_i - f(x_i, \theta))^2
$$

We simply "follow" the loss by deriving it with respect to the parameters $\theta$:

$$
\theta^{m+1} = \theta^m - \rho \nabla R(\theta^m)
$$

Where $\rho$ is the learning rate. And the gradient is computed as follows:

$$
\nabla R(\theta^m) = \frac{\partial R(\theta)}{\partial \theta} |_{\theta^m}
$$

This moves us towards the function's local minimum.

We might get stuck in a local minimum, but this is not a problem in practice.
There exist other optimization algorithms that can help us escape local minima.
In particular today, we will use the Adam optimizer which uses the moments of the gradients to escape the local minima.

### Learning rate

The gradient of the loss tells us the direction in which to take our next step.
It does not tell us how big the step should be.
The learning rate is a hyperparameter that needs to be tuned.
If it's too small, the algorithm will take too long to converge.
If it's too large, the algorithm might overshoot the minimum and not converge at all.

## Hyperparameters tuning

First of all we must choose an architecture: how many layers, nodes per layer, activation functions, etc.
Then we must choose the learning rate, the batch size, the number of epochs, etc.
The latter are called hyperparameters.

Hyperparameter have to the tuned in order to get the best performance out of the model.

This can be achieved by defining a search space and then using a search algorithm to find the best hyperparameters.
This can easily be achieved by employing a third-party library such as `optuna`.

There are two main techniques for hyperparameter tuning:

1. **Grid search**: we construct a grid of hyperparameters and evaluate the model for each combination. The model with the lowest validation loss is chosen.
    - **Drawbacks**: it is computationally expensive since we could have a huge grid.
2. **Random search**: we randomly chose some points in the search space and evaluate the model for each combination.
    The model with the lowest validation loss is chosen.
    - **Drawbacks**: it is computationally expensive since we could have a huge search space.

The random search is done iteratively.
We start from a very borad range and iteratively narrow it down to more interesting regions.
This can even be applied to the architecture of the model itself.

### Pruning

Pruning is a technique to optimize the tuning process.

Very plainly we 'prune', i.e. stop the training of a model if it is not performing well in the beginning.
This can save us a lot of time.
