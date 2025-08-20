# Research code repository

Here I collect and organise the code used to generate the results that end up in peer-reviewed journal articles.

Each subdirectory is associated one-to-one to a specific publication and it includes simulation scripts, postprocessing routines and plotting utilities for reproducibility and transparency.

The structure of the subdirectories is consistent throughout the repo and follows this scheme

```bash
short_title_of_pub/
├── doc/                # LaTeX source code and compiled PDF of the arXiv preprint 
├── inc/                # Include files collecting the source code in src/ 
├── sim/                # Scripts used to generate the figures in the paper 
└── src/                # Functions and modules used by the scripts in sim/ 
```

In each subdirectory a `README` file explains the main content of the paper, provides links to the open-access preprint and published version of record, shows some cool pictures/animations and provides instructions on how to cite the paper.

The subdirectories are not sorted in any meaningful order. 
Viewers interested in one particular work can follow one of the links listed below:
- [Localised patterns of the Lugiato-Lefever equation in 2 dimensions](./2d_localised_patterns_in_LLE/) _work in progress_;
- [Rate-induced tipping in a controlled model of game-theory](./R_tipping_in_game_theory/) _work in progress_;
- [Functional characterisation of rate-induced tipping](./rate_induced_tipping/) _work in progress_;
- [A probabilistic early-warning signal of saddle-node transitions from statistical mechanics](./statistical_ews_in_saddle_node_systems/) _work in progress_;
- My [PhD thesis](./phd_thesis/) _work in progress_.
