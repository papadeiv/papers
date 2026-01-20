# Localised patterns of the Lugiato-Lefever equation (LLE) in 2 dimensions

## 🚀 Description 

_Work in progress_

### ✏️📋 Outline

_Work in progress_

### 🔍📜 How to cite

The numerical method underlying the continuation of a discretised PDE is implemented by Daniele Avitabile. Please cite his work as

```bash
@misc{avitabile2020,
  author = {Daniele Avitabile},
  title = {{Numerical Computation of Coherent Structures in Spatially-Extended Systems}},
  month = May,
  year = 2020,
  publisher = {Zenodo},
  version = {v1.2-alpha},
  doi = {10.5281/zenodo.3821169},
  url = {https://doi.org/10.5281/zenodo.3821169}
}
```

## 📦 Structure of the repo

This repository is organised as follows

```bash
root/
├── inc/                # Include scripts importing local functions in src/
├── res/                # Figures and animations obtained by the results of sim/ 
├── sim/                # Actual numerical simulations of the work
└── src/                # Reusable functions implementing the algorithms used in sim/
```

### ⚙️🛠️ Organisation of the experiments

The actual numerical continuation algorithm for the pattern-forming PDE is implemented in MATLAB by [Daniele Avitabile](https://www.danieleavitabile.com/numerical-computation-of-coherent-structures-in-spatially-extended-systems/).

```bash
sim/
├── data/               # Results of the secant continuation of the files in main/ stored in the .mat format 
├── main                # Numerical continuation scripts and local copy of Daniele Avitabile's code
└── plotting            # Plotting routines of the continuation in 1 and 2 dimensions
```

### 💡🔬 What does each experiment do?

What follows is a brief description of the numerical experiments so that you can understand their purpose without interpreting it from the code:

- `fig:snaking_1d`: bifurcation diagram in the L2-norm of steady-state solutions of the LLE in 1 dimension showing homoclinic snaking [[Parra-Rivas et al. 2018]](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.97.042204). Candidate solutions (yellow stars) are selected to continue them in the second dimension via homotopy.

![](res/2.jpeg)

- `fig:homotopy_continuation_step_xxx`: continuation in 2 dimensions of selected solutions of the bifurcation diagram above

| Unsuccesfull homotopy | Succesfull homotopy |
| --------------------- | ------------------- |
| ![](res/3.jpeg)       | ![](res/4.jpeg)     |

- `fig:snaking_2d`: bifurcation diagram of the steady-state obtained via succesfull homotopy continuation from the 1-dimensional pattern (yellow star)

| Bifurcation diagram in 2d        |

| --------------- | -------------- |
| ![](res/5.jpeg) | ![](res/6.gif) |
