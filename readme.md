README
===============
**The repository is part of the publication: The choice of the objective function in flux balance analysis is crucial for predicting replicative lifespans in yeast, Schnitzer et al., 2022**

The model is adapted from our publication: Integrated model of yeast metabolism and replicative ageing reveals importance of metabolic phases during ageing, Schnitzer and Österberg et al., 2022, (https://github.com/cvijoviclab/IntegratedModelMetabolismAgeing)


**Abstract**

Flux balance analysis (FBA) models are a powerful tool to study genome-scale models of the cellular metabolism, based on finding the optimal flux distributions over the network. While the objective function is crucial for the outcome, its choice, even though motivated by evolutionary arguments, has not been directly connected to related measures. Here, we used an available multi-scale mathematical model of the cellular metabolism, damage accumulation and ageing, to systematically test the effect of commonly used objective functions on features of replicative ageing in budding yeast, such as the number of cell divisions and the corresponding time between divisions. The simulations confirmed that assuming maximal growth is essential for reaching realistic lifespans. 
The usage of the parsimonious solution or the additional maximisation of a growth-independent energy cost can improve lifespan predictions, explained by enhancing respiratory and antioxidative activity, using resources otherwise allocated to cellular growth, specifically in early life. 
Our work provides a new perspective on choosing the objective function in FBA models by connecting it to replicative ageing.


## Structure of the repository

- **Functions** contains all functions needed to run the model, as well as a schematic picture of the code structure. To run a lifespan simluation, only the file "IntegratedModel.jl" is needed to be imported.

- **ModelFiles** contains relevant model files for the modules of the integrated model: Boolean model for the signalling, FBA model for the central carbon metabolism.

- **Simulations** contains all simulations relevant in the publication. *FlexibilityGrid* contains simulations to generate Figure 1,4 and S1. *ParsimoniousCheck_enzymes* contains simulations to generate Figure 3C. *ParsimoniousCheck_fluxes* contains simulations to generate Figure 3B and D and S2. Lastly, *ParsimoniousCheck_normalisedFluxes* contains simulations to generate 2A and B, 3A, S3 and S4. In each folder, the file type defines the purpose of the file (executable or parameter definitions: .jl, result: .txt, plotting function .R, plots .svg). 

## Required software

The simulations were tested with the **Julia** programming language, version 1.6.1, including important packages MAT (reading in FBA model), JuMP (optimisation toolbox) and Gurobi (solver for linear programs).
