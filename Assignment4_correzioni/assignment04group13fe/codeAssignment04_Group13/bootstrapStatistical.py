import numpy as np


def bootstrapStatistical(numberOfSamplesToBootstrap, returns):
    """
    The function generates a bootstrap sample from the given returns.

    INPUTS:
        numberOfSamplesToBootstrap:  The number of samples to bootstrap
        returns:                     The returns from which to generate the bootstrap sample

    OUTPUTS:
        samples:  The bootstrap sample
    """
    indices = np.random.choice(np.arange(0, len(returns)-1), size=numberOfSamplesToBootstrap, replace=False)  # take the indices of the sample
    sorted_indices = np.sort(indices)  # order indices
    samples = returns.iloc[sorted_indices]  # take the returns to make the sample
    return samples

