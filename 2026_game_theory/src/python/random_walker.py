########################################################################################
#      Implementing a simple random walker under different sampling distribution       # 
########################################################################################

import numpy as np
from scipy import special 
from Process import Process 
from Estimator import Estimator

# Define different pdfs for the diffusion (residuals) of the stochastic process 
# Gaussian distribution
def gaussian(x:float):
    mu = 0.0
    sigma = 1
    return (1.0/(sigma*np.sqrt(2*np.pi)))*(np.exp(-0.5*np.square((x-mu)/sigma)))

# Gamma distribution
def gamma(x:float):
    k = 7.5
    theta = 1.0
    return (1.0/(special.gamma(k)*(np.power(theta, k))))*np.power(x, k-1.0)*np.exp(-x/theta)

# Beta distribution
def beta(x:float):
    a = 2.0
    b = 2.0
    return (1.0/special.beta(a, b))*np.power(x, a-1.0)*np.power(1.0-x, b-1.0)

# Define the deterministic drift (trend) of the stochastic process
def drift(x:float):
    return 1*x

# Define the periodic fluctuations (seasonality) of the stochastic process
def seasonality(x:float):
    return 1*((x/10)*np.sin(x)*np.cos(x) - np.sin(3*x))

# Create the stochastic processes from each of the sampling pdfs 
gaussian_rw = Process(density=gaussian, domain=(-10,10))
gamma_rw = Process(density=gamma, domain=(0,20))
beta_rw = Process(density=beta, domain=(0,1))

# Generate the timeseries of the processes with Nt realizations 
Nt = 1000
gaussian_rw.evolve(Nt, drift=drift, season=seasonality)

# Detrend the timeseries for better analysis
gaussian_rw.detrend(mode='EMD', order=1)

# Plot the histogram of the drawn samples
gaussian_rw.plot()

# Estimate and plot statistical indicators of the timeseries
ts = Estimator(gaussian_rw.detrended.ts, timescale=10)
ts.autocorrelation(10, lag=1)
