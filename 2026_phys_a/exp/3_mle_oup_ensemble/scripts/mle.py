from pymle.models import OrnsteinUhlenbeck
from pymle.sim.Simulator1D import Simulator1D
from pymle.fit.AnalyticalMLE import AnalyticalMLE
from pymle.core.TransitionDensity import EulerDensity 
import matplotlib.pyplot as plt
import numpy as np

def infer(sample, dt):
    # Define the regression model
    model = OrnsteinUhlenbeck()

    # Provide an initial guess and parameter space
    #guess = np.array([1, 0.1, 0.4])
    guess = np.random.rand(3)
    bounds = [(-10, 10), (-10, 10), (0.00001, 10)] 

    # Solve the MLE
    estimate = AnalyticalMLE(sample, bounds, dt, density = EulerDensity(model)).estimate_params(guess)

    return estimate.params 
