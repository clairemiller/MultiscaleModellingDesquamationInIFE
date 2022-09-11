# Multiscale modelling of desquamation in the interfollicular epidermis
Claire Miller, Edmund Crampin, James Osborne (2022). PLOS Computational Biology 18(8): e1010368. https://doi.org/10.1371/journal.pcbi.1010368

## Repository structure

This repository contains the code used for the results in the above manuscript. The repository folder structure is as follows:
- `src/` the extra class files required to run the simulation code using Chaste
- `test/` the test files specifying the simulation setup to run for the results in the paper
- `figures/` the R code to generate the figures for the paper using the results from the Chaste simulations.

## Running the simulations using Chaste

### Core Chaste code
This code should run using the core version of Chaste, which can be found [here](https://chaste.cs.ox.ac.uk/trac/wiki), along with documentation on how to install and run the code. 
In case there is any issue with backwards compatability in future versions of the code, the Chaste version used for the results in this manuscript can be pulled from [my Chaste repo](https://github.com/clairemiller/Chaste). 

### Simulations for the paper

The code must be compiled using the compile only flag `co=1`, and then simulations are run using the flag `-seed [seed]` to specify the seed for the simulation.
All results in the paper used seeds 0--9. 
I compiled the code using scons, an example of how to compile and run the simulations with different levels of inhibitor and a seed of 1 would be:
```
>>> cd [Chaste Directory]
>>> scons b=GccOpt co=1 projects/MultiscaleModellingDesquamationInIFE/test/TestDifferentLevelsInhibitor.hpp &> build_output.log
>>> ./projects/MultiscaleModellingDesquamationInIFE/build/optimised/TestDifferentLevelsInhibitorRunner -seed 1 &> run_output_seed01.log
```
Note, the simulations take multiple days.
Each test case is setup with multiple simulations to run consecutively (i.e. all 5 levels of inhibitor would be run with seed 1 consecutively using the above commands).
For each scenario, a fill simulation will be run first (code calls `RunFillTissue` in `test/FillTissueFunctionsProject3.hpp`), before the actual simulation is run to initialise the system with a filled tissue domain.

## Connecting the code with the paper figures
The following details which test and figures code to run to generate each of the results and figures found in the manuscript.

### Figures R code setup
The figures code assumes all results are contained in subdirectories of one main results directory. The path to the parent results directory needs to be specified in line 79 of `figures/scripts/paper_plot_theme.R`: `results_parent_dir <- "~/Chaste/Chaste_Results/MultiscaleModellingDesquamationInIFE/"`
