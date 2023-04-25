# Crekapp
This repository contains code to produce the results presented in the publication "Proteomics and constraint-based modelling reveal enzyme kinetic properties of Chlamydomonas reinhardtii on a genome scale"

## Dependencies
 - [R](https://www.r-project.org/) (tested on version 4.2.1)
 - [matlab](https://www.mathworks.com/help/install/install-products.html) (tested on version 2020b & 2022b)
 - [COBRA](https://opencobra.github.io/cobratoolbox/stable/installation.html) toolbox (tested on 2022 release)
 - [RAVEN] (https://github.com/SysBioChalmers/RAVEN) (tested on version 2.5.0)
 - [GECKO](https://github.com/SysBioChalmers/GECKO) toolbox (tested on version 2.0.2)
 - [Gurobi solver](https://www.gurobi.com/documentation/9.5/quickstart_mac/software_installation_guid.html) (tested on version 9.5.2)

## Setup (~ 30min)
After making sure all dependencies are installed and properly setup (see links above). Switch working directory to this git repository and open an R console and from within run
`source("Rsetup.r")` to install nescessary R packages.
Edit the `Matlab_startup.m` file to add the paths of the installed dependencies.

## Reproducing results (on linux machine - scripts produce verbose output to comman line no reason to worry)
1. Set this repository as working directory 

2. Run `Rscript Program/fit_chemostatdat.r` to obtain a model for maximum acetate uptake. (On windows it is sometimes only possible to source the scripts from within the R console using `source("Program/fit_chemostat.r")`)(~2-3s on Ryzen5 4000 16 GB RAM)

3. In R run `Rscript Program/QCC_smy_main2.r` to process the raw QCC data into and generate plots of QconCat proteomics overview statistics. (Also here on windows alternatively you can use `source("Program/QCC_smy_main2.r")`)(~90s on Ryzen5 4000 16 GB RAM)

4. In Matlab run `Matlab_startup` to set path for depenencies (edit as mentioned above) (1s on Ryzen5 4000 16 GB RAM)

5. In Matlab run `GECKO_startup()` to generate pcGEMs from autotrophic, mixotrohpic and heterotrophic chlamydomonas models. This will create a log file with GECKO output in the working directory. (~3h on  Ryzen5 4000 16 GB RAM)

6. In Matlab run `comp_ecModel_rescale` to generate metabolic model predictions from the GEM and pcGEMs. (~2m on Ryzen5 4000 16 GB RAM)

7. In R run `Program/plotprogrep_202207.r` to generate figures and statistics presented in the paper. (30s on Ryzen5 4000 16 GB RAM)
