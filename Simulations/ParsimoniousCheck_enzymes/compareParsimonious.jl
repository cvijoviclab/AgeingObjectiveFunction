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
outputFile = "enzymes.txt"
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
    shiftIdx = [maxTimesteps+1, maxTimesteps+1]

    # initialise enzymes and some more reations needed
    enzymes = cell.fba.model[:enzymes]
    nEnzymes = length(enzymes)
    reactionPathways = cell.fba.reactionPathways
    growthFlux = cell.relevantRefs["Growth"]
    glucoseUptakeFlux = cell.relevantRefs["Uptake of glucose"]
    ethanolUptakeFlux = cell.relevantRefs["Production of ethanol (reversible)"]
    ethanolProductionFlux = cell.relevantRefs["Production of ethanol"]

    # get the mapping between enzymes and reactions
    enzymesToComponents, componentsToFluxes = mapEnzymeIdx(cell.fba)

    # go through time steps individually and save enzymes over time
    # NOTE : since we go only one time step at a time the simulateLife!() 
    # function always start from 0 in each time step, while we count time and rls
    # and the phases here instead
    enzymeDynamics = zeros(Float64, nEnzymes + 1, maxTimesteps)
    for n = 1 : maxTimesteps

        # save the currently used inputs for the Boolean model
        tmpBooleanInputs = cell.latestBooleanInputs[1]

        # run one time step (maximalTimesteps = 1)
        status, tmpRls = simulateLife!(cell, mass, sizeProportion, retention, 
                                       timestep, 1, regulationFactor, objective1, 
                                       objective2)

        # break if there is no solution anymore
        !(status == "alive") ? break : nothing

        # save enzymes
        tmpEnzymes = value.(enzymes)
        enzymeDynamics[:, it] = [time; tmpEnzymes]

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
    # calculate the mean of all fluxes in the end of phase I and phase II
    # note that since we don't allow ethanol uptake there is no phase III here,
    # i.e. cell death means automatically end of phase II
    meanPhase1 = mean(enzymeDynamics[2:end, 1:shiftIdx[1]-1], dims = 2)
    stdPhase1 = std(enzymeDynamics[2:end, 1:shiftIdx[1]-1], dims = 2, 
                    mean = meanPhase1)
    meanPhase2 = mean(enzymeDynamics[2:end, shiftIdx[1]:it-1], dims = 2)
    stdPhase2 = std(enzymeDynamics[2:end, shiftIdx[1]:it-1], dims = 2, 
                    mean = meanPhase2)

    ############################################################################
    # prepare names of objectives
    if objective2.description != ""
        name1 = objective1.description * " (" * string(objective1.sense)[1:3] * ")"
        name2 = objective2.description * " (" * string(objective2.sense)[1:3] * ")"
    else
        name1 = objective1.description * " (" * string(objective1.sense)[1:3] * ")"
        name2 = ""
    end

    # add one row for each pathway each enzyme is involved in
    for e = 1:nEnzymes
        enzymeName = cell.fba.enzymeNames[e, 2]
        involvedReactions = componentsToFluxes[enzymesToComponents[e]]
        involvedPathways = unique(reactionPathways[involvedReactions])

        for pw in involvedPathways
            output[:, outputIdx] = [name1, name2, tmpParsimonious, "phase I",
                                    tmpFlex1, tmpFlex2, rls, averageGenerationTime,
                                    pw, enzymeName, meanPhase1[e], stdPhase1[e]]
            output[:, outputIdx + 1] = [name1, name2, tmpParsimonious, "phase II",
                                        tmpFlex1, tmpFlex2, rls, averageGenerationTime,
                                        pw, enzymeName, meanPhase2[e], stdPhase2[e]]
            outputIdx += 2
        end
    end
end

output = output[:, 1:outputIdx-1]


################################################################################
@printf("Save in output file ...\n")
################################################################################
pretext = 
"# ENZYME USAGES WITH OR WITHOUT PARSIMONIOUS SOLUTION
# FBA model $fbaPath
# Boolean model $booleanRulesPath, $booleanSpeciesPath
# TF layer $TFpath
# basal damage formation f0 = $formationRate0
# basal damage repair r0 = $repairRate0
objective 1\tobjective 2\tparsimonious\tphase\tflexibility 1\tflexibility 2\trls\t" *
"average generation time\tpathway\tenzyme\taverage usage\tstd usage
"

open(outputFile, "w") do file
    write(file, pretext)
    writedlm(file, permutedims(output), "\t")
end

