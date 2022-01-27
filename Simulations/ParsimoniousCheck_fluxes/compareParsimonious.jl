################################################################################
################################################################################
# simulate the lifespan of a cell with an integrated model 
# (ecFBA for the metabolism, Boolean for the signalling, ODE
# for damage accumulation) for maximal growth (with potential second optimisation)
# to find the difference between the parsimonious and non-parsimonious solution
################################################################################
################################################################################

using Printf
using DelimitedFiles
using Statistics
using StringDistances
include("../../Functions/IntegratedModel.jl")


################################################################################
@printf("Define parameters ... \n")
################################################################################
include("parameters.jl")

# we assume here that for each element in firstObjective there exists a corres-
# ponding element in firstSense, secondObjective and secondSense
parameterCombinationIdx = [[i, j, k, l] for i = 1:length(firstObjective), 
                                            j = 1:length(firstFlexibility),
                                            k = 1:length(secondFlexibility),
                                            l = 1:length(parsimonious)]
parameterCombinationIdx = parameterCombinationIdx[:]
nCombinations = size(parameterCombinationIdx, 1)


################################################################################
@printf("Initialise the output ... \n")
################################################################################
outputFile = "fluxes.txt"
output = Array{Any, 2}(undef, 12, nCombinations * 500)
outputIdx = 1


################################################################################
@printf("Go through parameters and simulate cell ... \n")
################################################################################
for i = 1 : nCombinations

    @printf("\tcombination %i/%i\n", i, nCombinations)

    global outputIdx

    ############################################################################
    # extract parameters for this simulation
    indices = parameterCombinationIdx[i]
    tmpFlex1 = firstFlexibility[indices[2]]
    tmpFlex2 = secondFlexibility[indices[3]]
    tmpParsimonious = parsimonious[indices[4]]

    # create the two objectives for this simulation
    objective1 = Objective(firstObjective[indices[1]], 
                           mapObjectiveToCoefficients(firstObjective[indices[1]]),
                           firstSense[indices[1]], tmpFlex1, tmpParsimonious)
    objective2 = Objective(secondObjective[indices[1]], 
                           mapObjectiveToCoefficients(secondObjective[indices[1]]),
                           secondSense[indices[1]], tmpFlex2, tmpParsimonious)

    ############################################################################
    # initialise cell
    cell = createCell(fbaPath, booleanSpeciesPath, booleanRulesPath, TFpath,
                      nDelay, maxGrowthRate, formationRate0, repairRate0,
                      P, D, mass, [glucoseThreshold, damageThreshold, 
                      trxThreshold], objective1, objective2)
    set_upper_bound(cell.relevantRefs["Production of ethanol (reversible)"],
                    0.0)
    simulateLife!(cell, mass, sizeProportion, retention, timestep, nDelay,
                  regulationFactor, objective1, objective2)
   
    # initialise time and division times
    time = 0.0
    divisionTimes = [0.0]
    rls = 0
    it = 1

    # prepare finding the phase shifts
    initialGrowth = maxGrowthRate
    phaseShift = [false, false]
    shiftIdx = [maxTimesteps + 1, maxTimesteps + 1]

    # go through time steps individually and save output fluxes over time
    # NOTE : since we go only one time step at a time the simulateLife!() 
    # function always start from 0 in each time step, while we count time and rls
    # and the phases here instead
    modelFluxes = cell.fba.model[:fluxes]
    growthFlux = cell.relevantRefs["Growth"]
    glucoseUptakeFlux = cell.relevantRefs["Uptake of glucose"]
    ethanolUptakeFlux = cell.relevantRefs["Production of ethanol (reversible)"]
    ethanolProductionFlux = cell.relevantRefs["Production of ethanol"]
    nFluxes = length(modelFluxes)
    fluxes = zeros(Float64, nFluxes + 1, maxTimesteps)
    for n = 1 : maxTimesteps

        # save the currently used inputs for the Boolean model
        tmpBooleanInputs = cell.latestBooleanInputs[1]

        # run one time step (maximalTimesteps = 1)
        status, tmpRls = simulateLife!(cell, mass, sizeProportion, retention, 
                                       timestep, 1, regulationFactor, objective1, 
                                       objective2)

        # break if there is no solution anymore
        !(status == "alive") ? break : nothing

        # save fluxes
        tmpFluxes = value.(modelFluxes)
        fluxes[:, it] = [time; tmpFluxes]

        # update rls and division times
        if tmpRls == 1
            rls += 1
            push!(divisionTimes, time)
        end

        # distinguish between phases (respect precision of solver)
        it == 1 ? initialGrowth = value(growthFlux) : nothing
        highestGrowth = value(growthFlux) - initialGrowth * 0.95 >= -precision
        ethanolUptake = (value(ethanolUptakeFlux) - value(ethanolProductionFlux)) > 
                        precision
        if !highestGrowth && !ethanolUptake && phaseShift[1] == false
            phaseShift[1] = true
            shiftIdx[1] = it
        elseif ethanolUptake && phaseShift[2] == false
            phaseShift[2] = true
            shiftIdx[2] = it
        end

        # update time
        time += timestep
        it += 1
    end  

    # get average generation time
    if rls > 0
        averageGenerationTime = mean(divisionTimes[2:end-1] .- divisionTimes[1:end-2])
    else
        averageGenerationTime = 0.0
    end

    ############################################################################
    # Go back to bi-directional network and remove isoenzymes
    reactionsToDelete = []

    # get all pseudo components that correspond to reactions arm reactions for 
    # isoenzymes or slack variables
    armComponents = findall(x -> contains(x, "pmet_") || contains(x, "pseudo"), 
                            cell.fba.componentNames)

    # go through the pseudo arm components and remove the reactions that are for the 
    # different enzymes (do not contain "arm" in their name)
    reactionNames = cell.fba.reactionNames
    reactionPathways = cell.fba.reactionPathways
    stochiometry = cell.fba.model[:stochiometry]
    for i = 1 : length(armComponents)
        component = armComponents[i]
        involvedFluxIdx = findall(x -> x != 0, 
                                  normalized_coefficient.(stochiometry[component], 
                                                          cell.fba.model[:fluxes]))
        involvedReactions = reactionNames[involvedFluxIdx]
        toDelete = .!contains.(involvedReactions, "(arm")
        reactionsToDelete = [reactionsToDelete; involvedFluxIdx[toDelete]]
    end
    reactionsToDelete = unique(reactionsToDelete)

    # remove the corresponding reactions from the output 
    reactionsToKeep = setdiff(1:nFluxes, reactionsToDelete)
    fluxes = fluxes[[1; reactionsToKeep.+1], :]
    reactionNames = reactionNames[reactionsToKeep]
    reactionPathways = reactionPathways[reactionsToKeep]

    # reset the vector for deletions for the reversible reactions
    reactionsToDelete = []

    # get all (remaining) reactions that are reversible
    reversibleReactions = findall(x -> contains(x, "(reversible)"), reactionNames)

    # go through the reversible reactions, find the forward reaction and use 
    # the sum of both reactions as the forward reaction
    for i = 1 : length(reversibleReactions)

        reversible = reactionNames[reversibleReactions[i]]

        # remove (reversible) in the string and find the closest reaction name
        # which corresponds to the forward reaction
        forward = replace(reversible, "(reversible)" => "")
        forwardReaction, idx = findnearest(forward, reactionNames, Levenshtein())

        # substract the flux of the reversible reaction from the forward reaction
        fluxes[idx+1, :] .-= fluxes[reversibleReactions[i]+1, :]
    end

    # delete the reversible reactions
    reactionsToKeep = setdiff(1:length(reactionNames), reversibleReactions)
    fluxes = fluxes[[1; reactionsToKeep.+1], :]
    reactionNames = reactionNames[reactionsToKeep]
    reactionPathways = reactionPathways[reactionsToKeep]

    ############################################################################
    # calculate the mean of all fluxes in the end of phase I and phase II
    # note that since we don't allow ethanol uptake there is no phase III here,
    # i.e. cell death means automatically end of phase II
    meanPhase1 = mean(fluxes[2:end, 1:shiftIdx[1]-1], dims = 2)
    stdPhase1 = std(fluxes[2:end, 1:shiftIdx[1]-1], dims = 2, mean = meanPhase1)
    meanPhase2 = mean(fluxes[2:end, shiftIdx[1]:it-1], dims = 2)
    stdPhase2 = std(fluxes[2:end, shiftIdx[1]:it-1], dims = 2, mean = meanPhase2)

    # save all remaining fluxes
    if objective2.description != ""
        name1 = objective1.description * " (" * string(objective1.sense)[1:3] * ")"
        name2 = objective2.description * " (" * string(objective2.sense)[1:3] * ")"
    else
        name1 = objective1.description * " (" * string(objective1.sense)[1:3] * ")"
        name2 = ""
    end
    nNewRows = length(reactionPathways)

    # for phase I
    outputParams = repeat([name1, name2, tmpParsimonious, "phase I", tmpFlex1, 
                           tmpFlex2, rls, averageGenerationTime], 1, nNewRows)
    outputFluxes = [reactionPathways reactionNames meanPhase1 stdPhase1] 
    output[:, outputIdx:outputIdx+nNewRows-1] = [outputParams; permutedims(outputFluxes)]
    outputIdx += nNewRows

    # for phase II
    outputParams = repeat([name1, name2, tmpParsimonious, "phase II", tmpFlex1, 
                           tmpFlex2, rls, averageGenerationTime], 1, nNewRows)
    outputFluxes = [reactionPathways reactionNames meanPhase2 stdPhase2] 
    output[:, outputIdx:outputIdx+nNewRows-1] = [outputParams; permutedims(outputFluxes)]
    outputIdx += nNewRows
end

output = output[:, 1:outputIdx-1]


################################################################################
@printf("Save in output file ...\n")
################################################################################
pretext = 
"# FLUXES WITH OR WITHOUT PARSIMONIOUS SOLUTION
# FBA model $fbaPath
# Boolean model $booleanRulesPath, $booleanSpeciesPath
# TF layer $TFpath
# basal damage formation f0 = $formationRate0
# basal damage repair r0 = $repairRate0
objective 1\tobjective 2\tparsimonious\tphase\tflexibility 1\tflexibility 2\trls\t" *
"average generation time\tpathway\treaction\taverage flux\tstd flux
"
open(outputFile, "w") do file
    write(file, pretext)
    writedlm(file, permutedims(output), "\t")
end

