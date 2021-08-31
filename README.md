# Increased Dynamic Range Example

This folder contains example code for applying a two-step reconstruction of
widely differing particle concentrations in magnetic particle imaging (MPI). 
The method combines two reconstruction with suitable parameters in order to cope with the different concentrations.

The method is described in the associated publication

M. Boberg, N. Gdaniec, P. Szwargulski, F. Werner, M. Möddel, and T. Knopp. (2021). Simultaneous imaging of widely differing particle concentrations in MPI: problem statement and algorithmic proposal for improvement. *Physics in Medicine & Biology*, 66(9), 095004. doi: [10.1088/1361-6560/abf202](https://iopscience.iop.org/article/10.1088/1361-6560/abf202).



## Installation

In order to use this code one first has to download [Julia](https://julialang.org/) (version 1.6 or later), clone this repository and navigate to the folder in the command line. The example script automatically activates the environment and installs all necessary packages.

## Execution
After installation the example code can be executed by running `julia` and entering
```julia
include("example.jl")
```
This will first download all data and then perform a reconstruction.
Parameters of the reconstruction are documented in the Julia script and can be
changed. After the reconstruction is done, the script will open a plotting window
and show the reconstruction result. 


## Open MPI Data

The measurement data associated to this project is about 2.1 GB large and will be downloaded and stored automatically, when the code is executed for the first time.
It is published under a [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/legalcode) license and can be found here:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5006979.svg)](https://doi.org/10.5281/zenodo.5006979)
