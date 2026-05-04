from pymle.models import OrnsteinUhlenbeck
from pymle.sim.Simulator1D import Simulator1D
from pymle.fit.AnalyticalMLE import AnalyticalMLE
from pymle.core.TransitionDensity import EulerDensity 
import matplotlib.pyplot as plt
import numpy as np

# ----- Parameters ----- #
X0 = 0.4
K = 3
μ = 0.3
σ = 0.2

model = OrnsteinUhlenbeck()
model.params = np.array([K, μ, σ])

# ----- Settings ------#
T = 5
ω = 250
dt = 1./ω
seed = 123

# ----- Simulation ------#
problem = Simulator1D(X0, T*ω, dt, model).set_seed(seed) 
sample = problem.sim_path()

# ----- MLE ----- #
Ω = [(0.01, 10), (0, 4), (0.01, 1)]
guess = np.array([1, 0.1, 0.4])
estimate = AnalyticalMLE(sample, Ω, dt, density = EulerDensity(model)).estimate_params(guess)
print(estimate)

# ----- Figure ------#
plt.plot(sample)
plt.show()
