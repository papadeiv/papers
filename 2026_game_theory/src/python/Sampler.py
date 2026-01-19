#######################################################################################################################################
#       This class implements a non-uniform random variate generator (https://en.wikipedia.org/wiki/Inverse_transform_sampling) 
#       to draw samples from a specific probability density function (pdf) provided as input by the user. The sampling method
#       that has been used is the inverse transformation method (https://en.wikipedia.org/wiki/Inverse_transform_sampling) which
#       is implemented by the stats module in SciPy and doesn't require the knowledge of the cumulative distribution (CDF)
#######################################################################################################################################

import numpy as np
from scipy import stats                                    # https://docs.scipy.org/doc/scipy/tutorial/stats.html
from scipy.stats import sampling                           # https://docs.scipy.org/doc/scipy/tutorial/stats/sampling.html 
from scipy.stats.sampling import NumericalInversePolynomial # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.sampling.NumericalInversePolynomial.html#scipy.stats.sampling.NumericalInversePolynomial 
from matplotlib import pyplot as plt                       # https://matplotlib.org/stable/users/indext
import plotly as ply                                       # https://plotly.com/python/
import seaborn as sb                                       # https://seaborn.pydata.org/tutorial.html

# Define the input pdf from which the target CDF is interpolated (https://en.wikipedia.org/wiki/List_of_probability_distributions)
class InputDensity:
    def __init__(self, function):
        self.density = function 

    def pdf(self, x:float):
        return self.density(x)

# Define the non-uniform generator
class Sampler:
    def __init__(self, function, domain=(float('-inf'),float('inf'))):
        # Initialise input pdf
        self.density = InputDensity(function)
        # Initialise the interpolated CDF
        self.sampler = NumericalInversePolynomial(self.density, domain=domain)
        # Initialise the sampling seed
        self.seed = 1

    def sample(self, ndraws):
        # Generate uniformly distributed samples in [0,1]
        seed = np.random.default_rng(seed=self.seed)
        self.samples = np.empty([ndraws,1])
        for _ in range(ndraws):
            self.samples[_] = seed.random()

        # Draw samples from the interpolated CDF
        self.sampler.set_random_state(self.seed)
        self.draws = self.sampler.rvs(ndraws)
        # Update the sampling seed for the next run
        self.seed += 1
        return

    def plot(self):
        # Create a figure
        fig = plt.figure(figsize=[12.8,9.6], dpi=200, layout='tight')
        # Plot the histogram of the original samples drawn from U([0,1])
        ax1 = plt.subplot2grid((4,4), (0,0), colspan=1, rowspan=2, fig=fig)
        ax1.hist(self.samples, bins=10, color="orange", edgecolor="black", orientation='horizontal')
        # Plot the histogram of the samples drawn from the preimages of the CDF
        ax2 = plt.subplot2grid((4,4), (2,1), colspan=3, rowspan=2, fig=fig)
        ax2.hist(self.draws, bins=50, density=True, edgecolor="black")
        # Compute and plot the CDF from the samples drawn from the inversion method
        x = np.sort(self.draws)
        it, cumsum = 0,0
        cdf = np.zeros(self.draws.shape)
        for sample in x:
            if it < (self.draws.shape)[0] -1:
                dx = np.abs(sample - x[it+1])
            else:
                dx = 1/((self.draws.shape)[0] - 1)

            cumsum += self.density.pdf(sample)*dx
            cdf[it] = cumsum
            it += 1

        ax3 = plt.subplot2grid((4,4), (0,1), colspan=3, rowspan=2, sharex=ax2, fig=fig)
        ax3.plot(x, cdf, color = 'black')
        # Show the figure
        plt.show()
        return
