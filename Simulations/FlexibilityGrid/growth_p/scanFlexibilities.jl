################################################################################
################################################################################
# simulate the lifespan of a cell with an integrated model 
# (ecFBA for the metabolism, Boolean for the signalling, ODE
# for damage accumulation) for a grid of different objectives, flexibilities
# and regulation
################################################################################
################################################################################

using Printf
using DelimitedFiles
using Statistics
include("../../../Functions/IntegratedModel.jl")


################################################################################
@printf("Define parameters ... \n")
################################################################################
include("parameters.jl")

# we assume here that for each element in firstObjective there exists a corres-
# ponding element in firstSense, secondObjective and secondSense
parameterCombinationIdx = [[i, j, k, l] for i = 1:length(firstObjective), 
                                            j = 1:length(regulationFactor),
                                            k = 1:length(firstFlexibility),
                                            l = 1:length(secondFlexibility)]
parameterCombinationIdx = parameterCombinationIdx[:]
nCombinations = size(parameterCombinationIdx, 1)


################################################################################
@printf("Initialise the output ... \n")
################################################################################
outputFile = "flexibilityGrid.txt"
output = Array{Any, 2}(undef, 16, nCombinations)


################################################################################
@printf("Go through flexibility parameters and simulate cell ... \n")
################################################################################
for i = 1 : nCombinations

    @printf("\tcombination %i/%i\n", i, nCombinations)

    # extract parameters for this simulation
    indices = parameterCombinationIdx[i]
    tmpRegulation = regulationFactor[indices[2]]
    tmpFlex1 = firstFlexibility[indices[3]]
    tmpFlex2 = secondFlexibility[indices[4]]

    # create the two objectives for this simulation
    objective1 = Objective(firstObjective[indices[1]], 
                           mapObjectiveToCoefficients(firstObjective[indices[1]]),
                           firstSense[indices[1]], tmpFlex1, parsimonious)
    objective2 = Objective(secondObjective[indices[1]], 
                           mapObjectiveToCoefficients(secondObjective[indices[1]]),
                           secondSense[indices[1]], tmpFlex2, parsimonious)

    # initialise cell 
    cell = createCell(fbaPath, booleanSpeciesPath, booleanRulesPath, TFpath,
                      nDelay, maxGrowthRate, formationRate0, repairRate0,
                      P, D, mass, [glucoseThreshold, damageThreshold, 
                      trxThreshold], objective1, objective2)

    # do not allow ethanol uptake
    set_upper_bound(cell.relevantRefs["Production of ethanol (reversible)"],
                    0.0)

    # iterate nDelay times to fill up the recent Boolean input states
    simulateLife!(cell, mass, sizeProportion, retention, timestep, nDelay,
                  tmpRegulation, objective1, objective2)

    # simulate life
    results = simulateLife!(cell, mass, sizeProportion, retention, timestep, 
                            maxTimesteps, tmpRegulation, objective1, objective2)
    status, rls, averageGenerationTime, phaseSplit = results

    # save in output
    if objective2.description != ""
        name = objective1.description * " (" * string(objective1.sense)[1:3] * ") + " *
               objective2.description * " (" * string(objective2.sense)[1:3] * ")"
    else
        name = objective1.description * " (" * string(objective1.sense)[1:3] * ")"
    end

    output[:, i] = [name; tmpFlex1; tmpFlex2; tmpRegulation; rls; 
                    averageGenerationTime; phaseSplit[1];
                    phaseSplit[2]; phaseSplit[3]; status]
end


################################################################################
@printf("Save in output file ...\n")
################################################################################
pretext = 
"# FLEXIBILITY PARAMETER SCAN
# FBA model $fbaPath
# parsimonious solution $parsimonious
# Boolean model $booleanRulesPath, $booleanSpeciesPath
# TF layer $TFpath
# basal damage formation f0 = $formationRate0
# basal damage repair r0 = $repairRate0
objective\tflexibility 1\tflexibility 2\tregulation factor\trls\t" *
"average generation time\ttime phase 1\trls phase 1\tdamage phase 1\t"*
"time phase 2\trls phase 2\tdamage phase 2\ttime phase 3\trls phase 3\t" *
"damage phase 3\tstatus
"
open(outputFile, "w") do file
    write(file, pretext)
    writedlm(file, permutedims(output), "\t")
end

