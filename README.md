# Multiscale modelling of desquamation in the interfollicular epidermis
Claire Miller, Edmund Crampin, James Osborne (2022). PLOS Computational Biology 18(8): e1010368. https://doi.org/10.1371/journal.pcbi.1010368

## Repository structure

This repository contains the code used for the results in the above manuscript. The repository folder structure is as follows:

- `src/` the extra class files required to run the simulation code using Chaste
- `test/` the test files specifying the simulation setup to run for the results in the paper
- `figures/` the R code to generate the figures for the paper using the results from the Chaste simulations.
- `matlab_code/` the matlab code for the single cell solutions.

## Running the simulations using Chaste

### Core Chaste code
This code should run using the core version of Chaste, which can be found [here](https://chaste.cs.ox.ac.uk/trac/wiki), along with documentation on how to install and run the code. 
Clone this repository into the `projects` subdirectory of the Chaste code. 
In case there is any issue with backwards compatability in future versions of the code, the Chaste version used for the results in this manuscript can be pulled from [my Chaste repo](https://github.com/clairemiller/Chaste), commit `e254861c657d325872aba156cc39a8ae56210a4b`. 

### Simulations for the paper

The code must be compiled using the compile only flag `co=1`, and then simulations are run using the flag `-seed [seed]` to specify the seed for the simulation.
All results in the paper used seeds 0--9. 
I compiled the code using scons, an example of how to compile and run the simulations with different levels of inhibitor and a seed of 1 would be:

```
~: cd [Chaste Directory]
~/Chaste: scons b=GccOpt co=1 projects/MultiscaleModellingDesquamationInIFE/test/TestDifferentLevelsInhibitor.hpp &> build_output.log
~/Chaste: ./projects/MultiscaleModellingDesquamationInIFE/build/optimised/TestDifferentLevelsInhibitorRunner -seed 1 &> run_output_seed01.log
```
Note, the simulations take multiple days.
Each test case is setup with multiple simulations to run consecutively (i.e. all 5 levels of inhibitor would be run with seed 1 consecutively using the above commands).
For each scenario, a fill simulation will be run first (code calls `RunFillTissue` in `test/FillTissueFunctionsProject3.hpp`), before the actual simulation is run to initialise the system with a filled tissue domain.

## Generating the paper figures
The following details which test and figures code to run to generate each of the results and figures found in the manuscript.

### R code preliminary setup
The figures code assumes all results are contained in subdirectories of one main results directory. 
To generate the figures, first do the following:

1. The path to the parent results directory needs to be specified in line 79 of `figures/scripts/paper_plot_theme.R`: `results_parent_dir <- "~/Chaste/Chaste_Results/MultiscaleModellingDesquamationInIFE/"`
2. Open the R project `figures/MultiscaleModellingDesquamationInIFE.Rproj` in RStudio if using
3. If desired, modify the plotting theme functions and colour sets in `scripts/paper_plot_theme.R`

### Figure 2: single cell results
These results were simulated using the Matlab script `matlab_code/solve_movingcell.m` which generates a csv file: `matlab_code/data/movingcell_solutions.csv`.
To generate the plots from this code, run `figures/scripts/moving_cell_solutions.R`, modifying line 4 (`filename <- [path]/matlab_code/data/movingcell_solutions.csv`) for the correct csv output path.
The figure 2 plot file `figures/movingcell_effective.pdf` will be generated. 
A second plot file `figures/movingcell_actual.pdf` will also be generated, which is Figure E in S1 Text.

### Figure 4: normal homeostatic tissue
These results are for simulations where we assume the tissue is in normal homeostasis.
Results can be used from either simulation code `test/TestDifferentLevelsInhibitor.hpp` with 100% inhibitor levels, or `test/TestMutations.hpp` with 0 diseased cells.
Results will be the same using either set.
To generate the plots in this figure, run `figures/scripts/normal_system.R`.
The four plot figure files (Figure 4B--E) will be generated as `figures/healthy_[height/reactants/velocity/turnovertime].pdf`.

### Figures 5 and 6: different levels of inhibitor
These results are for simulations where we change the level of inhibitor each cell secretes.
Results come from the simulation output of `test/TestDifferentLevelsInhibitor.hpp`.
To generate the plots in these figures run `figures/scripts/vary_inhibitor.R`.
The six plot figure files will be generated as `figures/varyinhibitor_[desciptor].pdf`.

### Figure 7: heterogeneous recovery
These results are for the simulations where we set different proportions of cells to have no or full levels of inhibitor.
Results come from the simulation output of `test/TestMutations.hpp`.
To generate the plots in these figures run `figures/scripts/treated_system.R`.
The four plot figure files will be generated as `figures/treated_[descriptor].pdf`. 

### Figure 8: force plot
This is a plot of the force functions used in the simulations.
It shows the force as a function of cell separation distance, as calculated in `src/Forces/PalssonAdhesionForce.cpp`.
To regenerate this figure, run the code in `figures/scripts/forces_diagram.R`.
The figure will appear in `figures/force.pdf`.
