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
    guess = np.array([5, 0.5, 0.5])
    bounds = [(0, 12), (-1, 1), (0, 10)] 

    # Solve the MLE
    estimate = AnalyticalMLE(sample, bounds, dt, density = EulerDensity(model)).estimate_params(guess)

    return estimate.params 
