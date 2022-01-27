README
===============
**As part of publication: The choice of the objective function in flux balance analysis models is crucial for predicting replicative lifespans in yeast, Schnitzer et al., 2022**

The model is adapted from the publication: Integrated model of yeast metabolism and replicative ageing reveals importance of metabolic phases during ageing, Schnitzer and Österberg et al., 2022, (https://github.com/cvijoviclab/IntegratedModelMetabolismAgeing)

**Abstract**




## Structure of the repository

- **Functions** contains all functions needed to run the model, as well as a schematic picture of the code structure. To run a lifespan simluation, only the file "IntegratedModel.jl" is needed to be imported.

- **ModelFiles** contains relevant model files for the modules of the integrated model: Boolean model for the signalling, FBA model for the central carbon metabolism.

- **Simulations** contains all simulations relevant to the publication. *FlexibilityGrid* contains simulations to generate Figure 1,4 and S1. *ParsimoniousCheck_enzymes* contains simulations to generate Figure 3C. *ParsimoniousCheck_fluxes* contains simulations to generate Figure 3B and D and S2. Lastly, *ParsimoniousCheck_normalisedFluxes* contains simulations to generate 2A and B, 3A, S3 and S4. In each folder, the file type defines the purpose of the file (executable or parameter definitions: .jl, result: .txt, plotting function .R, plots .svg). 

## Required software

The simulations were tested with the **Julia** programming language, version 1.6.1, including important packages MAT (reading in FBA model), JuMP (optimisation toolbox) and Gurobi (solver for linear programs).
