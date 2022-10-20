# Crekapp
This repository contains code to produce the results presented in the publication "Proteomics and constraint-based modelling reveal enzyme kinetic properties of C. reinhardtii on a genome-scale"

## Dependencies
 - [R](https://www.r-project.org/) (tested on version 4.2.1)
 - [matlab](https://www.mathworks.com/help/install/install-products.html) (tested on version 2020b & 2022b)
 - [COBRA](https://opencobra.github.io/cobratoolbox/stable/installation.html) toolbox (tested on 2022 release)
 - [GECKO](https://github.com/SysBioChalmers/GECKO) toolbox (tested on version 2.0.2)
 - [Gurobi solver](https://www.gurobi.com/documentation/9.5/quickstart_mac/software_installation_guid.html) (tested on version 9.5.2)

## Setup

After making sure all dependencies are installed and properly setup (see links above). Switch working directory to this git repository and run `Rsetup.r` to install nescessary R packages.
Edit the `Matlab_startup.m` file to add the paths of the installed dependencies.

## Reproducing results (on linux machine)
1. Set this repository as working directory 

2. Run `Program\fit_chemostatdat.r` to obtain a model for maximum acetate uptake.

3. Run `Program\Program/QCC_smy_main2.r` to process the raw QCC data into and generate plots of QconCat proteomics overview statistics.

4. 