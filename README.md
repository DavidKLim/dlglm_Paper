# dlglm_Paper
This repository contains the code to reproduce all of the analyses, figures, and tables in the paper ``Deeply-Learned Generalized Linear Models with Missing Data".

Python module dependencies:

torch=1.3.1; torchvision=0.4.2; numpy=1.17.4; tensorflow = 1.14.0; tensorflow-probability=0.7.0; tensorboard=1.14.0; cloudpickle=1.2.2; pandas=0.25.2

## Installing dlglm R package
First, install the R package dlglm. Note that you may have to double check your current installation of Python, Pytorch, CUDA, and relevant Python modules before proceeding. Please reference the README in the package repository for more details.
```
devtools::install_github("https://www.github.com/DavidKLim/dlglm")
```

## Running dlglm analyses
Analyses on simulated data, UCI real data with simulated missingness, and the Bank Marketing data can be replicated by running the R scripts `runComparisons_sim.R`, `runComparisons_UCI.R`, and `runComparisons_BANKMARKETING.R`, respectively. Please note that these analyses, especially the simulation analyses, are computationally intensive. It is therefore highly recommended to run these analyses on a high performance computing cluster.


## Reproducing results
The figures and tables pertaining to the analyses can be reproduced by running the `processResults_sim.R`, `processResults_UCI.R`, and `processResults_BANKMARKETING.R` scripts. One supplemental figure can also be reproduced by running the `plot_cor_drybean.R` script.
