# Characterisation of rate-induced tipping on parameteric families of parameter shifts

## 🚀 Description 

Rate-induced tipping (R-tipping) is a mechanism of critical transitions in non-autonomous dynamical systems [[Ashwin et al. 2012]](https://royalsocietypublishing.org/rsta/article-abstract/370/1962/1166/114607/Tipping-points-in-open-systems-bifurcation-noise?redirectedFrom=fulltext).
Contrary to bifurcation-induced (B-tipping), the critical transition is manifested by the pullback attractor not being able to track the drift of a quasi-steady equilibrium (QSE) escaping its basing of attraction [[Ashwin et al. 2017]](https://iopscience.iop.org/article/10.1088/1361-6544/aa675b).
The drift of the QSE in _phase space_ is associated to the rate at which the parameter of a dynamical system is changes.
As such, the evolution law (__shift__ of the paramater) could reach a critical rate that triggers R-tipping in these non-autonomous systems.

### 🔬 Outline

In this project we attempt to numerically investigate how the functional properties of non-monotonic, parametric shifts, characterise R-tipping in 1-dimensional, non-autonomous, nonlinear dynamical systems.

This project is currently at a _preliminary investigation_ stage!

## 📦 Structure of the repo

This repository is organised as follows

```bash
root/
├── doc/                # LaTeX source code and compiled PDF of the preliminary report 
├── inc/                # Includes scripts importing local functions in src/
├── sim/                # Actual numerical simulations of the work
└── src/                # Reusable functions implementing the algorithms used in sim/
```

### 🛠️ Organisation of the experiments

There are many numerical experiments in this work. Each of them is made of a `main` script, a `postprocessing` script and a `plotting` script. Therefore the `sim` subdirectory is structured as follows

```bash
sim/
├── main/               # Computes solutions of the system under investigation and exports the resulting data (users should execute these scripts) 
│   ├── exp_1.jl         
│   ├── exp_2.jl        
│   ├── ...        
│   └── exp_n.jl      
├── plotting            # Plots and exports the figures relevant to the analysed data in postprocessing 
│   ├── exp_1.jl       
│   ├── exp_2.jl      
│   ├── ...        
│   └── exp_n.jl
└── postprocessing      # Postprocess the data exported by the main script and generates the analysis 
    ├── exp_1.jl    
    ├── exp_2.jl   
    ├── ...        
    └── exp_n.jl
```

### 💡 What does each experiment do?

What follows is a brief description of the numerical experiments so that you can understand their purpose without interpreting it from the code:

- `figure_01`: given the non-autonomous system in [[Ashwin et al. 2017]](https://iopscience.iop.org/article/10.1088/1361-6544/aa675b), we perform a sweep on the Lipschitz constant of the parameter shift to identify the region for which rate-induced tipping occurs.

![](doc/figures/rate_sweep.png)

- `figure_02` to `figure_06`: we propagate the non-autonomous system above, forward in time but with different families of parameter shifts (mostly being polynomials in tanh(t)) to attempt to characterise rate-induce tipping.

| Sim 2                      | Sim 3                      | Sim 4.a                      |
| -------------------------- | -------------------------- | ---------------------------- |
| ![](doc/figures/sim_2.png) | ![](doc/figures/sim_3.png) | ![](doc/figures/sim_4.a.png) |

| Sim 6.b                      | Sim 6.c                      | Sim 6.d                      |
| ---------------------------- | ---------------------------- | ---------------------------- |
| ![](doc/figures/sim_6.b.png) | ![](doc/figures/sim_6.c.png) | ![](doc/figures/sim_6.d.png) |
