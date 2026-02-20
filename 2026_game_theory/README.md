# Rate-induced tipping in a feedback-controlled, game-theoretic model

## 🚀 Description 

_Work in progress_

### ✏️  Outline

_Work in progress_

### 📜 How to cite

_Work in progress_

## 📦 Structure of the repo

This repository is organised as follows

```bash
root/
├── doc/                # LaTeX source code and compiled PDF of the arXiv preprint
├── inc/                # Include scripts importing local functions in src/
├── sim/                # Actual numerical simulations of the work
└── src/                # Reusable functions implementing the algorithms used in sim/
```

### ⚙️  Organisation of the experiments

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

- `fig:autonomous`: we want to simulate and characterise the different behaviour of the same initial conditions (ICs) of the _autonomous_ replicator equation proposed in [[Zino et al. 2025]](https://ieeexplore.ieee.org/document/10988641) in 4 different regimes separated by transcritical bifurcations in the 2-dimensional parameter space M = (α, β) ⊆ ℝ2.

| Solutions in different regimes                | Bifurcation set                            |
| ----------------------------------------------| -------------------------------------------|
| ![](doc/figures/fig:autonomous_solutions.png) | ![](doc/figures/fig:autonomous_bifset.png) |

- `fig:nonautonomous`: we introduce a time-changing parameter λ according to a shift law Λ(t) that makes the system non-autonomous. We investigate the tipping behavious of the solutions as the parameter is shifted in M however ruling out rate-induced tipping.

| Tipping solutions (red curves) at different rates (colorbar) of the parameter shift (dashed black curve)                                                            | | | |
| ---------------------------------------- | ---------------------------------------- | ---------------------------------------- | ---------------------------------------- |
| ![](doc/figures/fig:nonautonomous_1.png) | ![](doc/figures/fig:nonautonomous_2.png) | ![](doc/figures/fig:nonautonomous_3.png) | ![](doc/figures/fig:nonautonomous_4.png) |
