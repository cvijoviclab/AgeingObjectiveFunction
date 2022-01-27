################################################################################
################################################################################
# InitialiseJuMPModel.jl includes functions for initialising a JuMP model
# from a enzyme-constraint FBA model:
#   - outputModel = initialiseJuMPModel(filepath)
#   - prepareChemostatExperiment!(fba)
#   - prepareLifespanExperiment!(fba, maxGrowthRate)
#   - enzymeCombinations = getIsoenzymes(fba)
#   - coefficients = mapObjectiveToCoefficients(objective:)
################################################################################
################################################################################


# set Gurobi parameters
const precision = 0.000000001
const outputFlag = 0
const feasibilityTol = precision
const optimalityTol = precision
const mipGap = 0.000000000001
const presolve = 2
const timeLimit = 1000
const method = -1

# set a scaling factor for the enzyme constraints in order to avoid numerical 
# issues
const scalingFactor = 1000


################################################################################
# define the model struct
################################################################################
struct FBAModel
    model::Model
    reactionNames::Array{String, 1}
    componentNames::Array{String, 1}
    enzymeNames::Array{String, 2}
    reactionPathways::Array{String, 1}
    enzymePathways::Array{String, 1}
end


################################################################################
# FBA model initialisation
#
# input parameters:
#   - filepath : string with the path to the model file (.mat)
#
# output parameters:
#   - outputModel : FBA struct with the JuMP optimisation problem and all names
#                   of the model's reactions and components
# 
################################################################################
function initialiseJuMPModel(filepath::String)::FBAModel

    # read file
    file = matread(filepath)
    modelName = collect(keys(file))[1]
    model = file[modelName]
    
    # extract relevant model information
    componentNames = string.(model["metNames"][:])
    reactionNames = string.(model["rxnNames"][:])
    reactionPathways = string.(model["rxnPathways"][:])
    stochiometry = model["S"] 
    changes = model["b"]
    nComponents, nReactions = size(stochiometry)
    lowerBounds = model["lb"][:]
    upperBounds = model["ub"][:]
    enzymeNames = string.(model["enzymes"])
    enzymePathways = string.(model["enzymePathways"][:])

    # divide the reactions in actual reaction fluxes, enzyme mass pseudo
    # reactions and a protein pool pseudo reaction
    nEnzymes = sum(contains.(reactionNames, "draw_prot"))
    nFluxes = nReactions - nEnzymes - 1
    
    # set the solver environment and create a JuMP model
    optModel = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(),
                     "OutputFlag" => outputFlag,
                     "FeasibilityTol" => feasibilityTol,
                     "OptimalityTol" => optimalityTol,
                     "MIPGap" => mipGap,
                     "Presolve" => presolve,
                     "TimeLimit" => timeLimit,
                     "Method" => method))
    
    # initilise the variable fluxes with their corresponding bounds
    # here we assume that the reactions are ordered according to 
    # first fluxes, then enzymes and last the enzyme pool
    @variable(optModel, lowerBounds[i] <= fluxes[i = 1:nFluxes] <= upperBounds[i])
    @variable(optModel, lowerBounds[i+nFluxes] <= enzymes[i = 1:nEnzymes] <=
                        upperBounds[i+nFluxes])
    @variable(optModel, lowerBounds[end] <= pool <= upperBounds[end])
    
    # rescale the enzymatic constraints to avoid numerical issues
    enzymeComponents = contains.(componentNames, "prot_") .& 
                       .!contains.(componentNames, "_pool")
    stochiometry[enzymeComponents, :] .*= scalingFactor
    changes[enzymeComponents] .*= scalingFactor

    # add constraints from stochiometry
    @constraint(optModel, stochiometry,
                stochiometry * [fluxes; enzymes; pool] .== changes)

    # return
    return FBAModel(optModel, reactionNames, componentNames, enzymeNames, 
                    reactionPathways, enzymePathways)
end


################################################################################
# set manual constraints for the chemostat experiment (switch between 
# respiration and fermentation)
#
# input parameters:
#   - fba : FBA model struct with all important information
#
# according to Österberg et al., PLOS Computational Biology 2021
# with extended constraints on damage parts
################################################################################
function prepareChemostatExperiment!(fba::FBAModel)::Nothing
 
    # extract needed model features
    componentNames = fba.componentNames
    reactionNames = fba.reactionNames
    fluxes = fba.model[:fluxes]
    pool = fba.model[:pool]
    stochiometry = fba.model[:stochiometry]

    # extract relevant reactions
    acetateProduction = findfirst(x -> x == "Production of acetate", 
                                  reactionNames)
    pyruvateProduction = findfirst(x -> x == "Production of pyruvate", 
                                   reactionNames)
    glycerolProduction = findfirst(x -> x == "Production of glycerol", 
                                   reactionNames)
    
    # update the upper bounds
    curatedIndices = [acetateProduction, pyruvateProduction, glycerolProduction]
    curatedUpperBounds = [0.682, 0.055, 0.165]
    set_upper_bound.(fluxes[curatedIndices], curatedUpperBounds)

    # update also the protein pool
    # set to 0.038 since it fits the best for the chemostat experiment
    # (independent of the added damage reactions including extra enzymes: 
    #   without 0.46 * 0.1714 * 0.481964 ≈ 0.038
    #   with 0.46 * 0.1799 * 0.4592 ≈ 0.038)
    set_upper_bound(pool, 0.46 * 0.1714 * 0.481964)
    
    # set initial value for GAM
    growthReaction = findfirst(x -> x == "Growth", reactionNames)
    GAMComponent = findfirst(x -> x == "Maintainance for growth [cytosol]", 
                             componentNames)
    set_normalized_coefficient(stochiometry[GAMComponent],
                               fluxes[growthReaction], -18)

    # uncontrain glucose uptake
    glucoseUptake = findfirst(x -> x == "Uptake of glucose", reactionNames)
    set_lower_bound(fluxes[glucoseUptake], 0.0)
    set_upper_bound(fluxes[glucoseUptake], 1000.0)

    # constrain direct damage production (to avoid that all fluxes go through
    # those reactions since they are non-enzymatic)
    directDamageReactions = ["hydroxyl production via peroxynitrite (cytosol)", 
                             "hydroxyl production via peroxynitrite (mitochondria)"]
    directDamageFluxes = fluxes[map(x -> findfirst(y -> y == x, reactionNames), 
                                directDamageReactions)]
    set_upper_bound.(directDamageFluxes, 0.1) 

    # restrict that reaction to a small flux (needed to allow at least some flux
    # through the cytosolic Fenton-Haber-Weiss reactions)
    superoxideProduction = findfirst(x -> x == "superoxide production (cytosol)",
                                     reactionNames)
    set_upper_bound(fluxes[superoxideProduction], 0.001)

    # set an upper bound to the non growth associated maintenance value
    # to assume a young cell with low damage levels
    NGAM = fluxes[findfirst(x -> x == "NGAM", reactionNames)]
    set_lower_bound(NGAM, 0.001)
    set_upper_bound(NGAM, 0.01)

    # return
    return nothing
end


################################################################################
# set manual constraints for the lifespan experiment
#
# input parameters:
#   - fba : FBA model struct with all important information
#   - maxGrowthRate : upper bound on the growth rate 
#
################################################################################
function prepareLifespanExperiment!(fba::FBAModel, 
                                    maxGrowthRate::Float64)::Nothing
 
    # extract needed model features
    componentNames = fba.componentNames
    reactionNames = fba.reactionNames
    fluxes = fba.model[:fluxes]
    pool = fba.model[:pool]
    stochiometry = fba.model[:stochiometry]

    # update further upper bounds (taken from Leupold et al., 2019)
    acetateProduction = findfirst(x -> x == "Production of acetate", 
                                  reactionNames)
    pyruvateProduction = findfirst(x -> x == "Production of pyruvate", 
                                   reactionNames)
    glycerolProduction = findfirst(x -> x == "Production of glycerol", 
                                   reactionNames)
    curatedIndices = [acetateProduction, pyruvateProduction, glycerolProduction]
    curatedUpperBounds = [1.25, 0.1, 1.2]
    set_upper_bound.(fluxes[curatedIndices], curatedUpperBounds)

    # limit growth to maximal growth (measured by Linnea)
    growthReaction = findfirst(x -> x == "Growth", reactionNames)
    set_lower_bound(fluxes[growthReaction], 0.0)
    set_upper_bound(fluxes[growthReaction], maxGrowthRate)

    # update the protein pool
    # set to 0.038 since it fits the best for the chemostat experiment
    # (independent of the added damage reactions including extra enzymes: 
    #   without 0.46 * 0.1714 * 0.481964 ≈ 0.038
    #   with 0.46 * 0.1799 * 0.4592 ≈ 0.038)
    set_upper_bound(pool, 0.46 * 0.1799 * 0.4592)

    # unconstrain glucose uptake
    glucoseUptake = findfirst(x -> x == "Uptake of glucose", reactionNames)
    set_lower_bound(fluxes[glucoseUptake], 0.0)
    set_upper_bound(fluxes[glucoseUptake], 1000.0)

    # allow ethanol influx (motivated by Leupold et al., 2019)
    ethanolInflux = findfirst(x -> x == "Production of ethanol (reversible)", 
                              reactionNames)
    set_upper_bound(fluxes[ethanolInflux], 1000.0)

    # constrain direct damage production (to avoid that all fluxes go through
    # those reactions since they are non-enzymatic)
    directDamageReactions = ["hydroxyl production via peroxynitrite (cytosol)", 
                             "hydroxyl production via peroxynitrite (mitochondria)"]
    directDamageFluxes = fluxes[map(x -> findfirst(y -> y == x, reactionNames), 
                                directDamageReactions)]
    set_upper_bound.(directDamageFluxes, 0.1) 

    # restrict that reaction to a small flux (needed to allow at least some flux
    # through the cytosolic Fenton-Haber-Weiss reactions)
    superoxideProduction = findfirst(x -> x == "superoxide production (cytosol)",
                                     reactionNames)
    set_upper_bound(fluxes[superoxideProduction], 0.001)

    # set maximal value for NGAM (value from ecYeast 8, Lu et al., 
    # Nature Communications 2019)
    NGAM = fluxes[findfirst(x -> x == "NGAM", reactionNames)]
    set_upper_bound(NGAM, 0.7)

    # return
    return nothing
end


################################################################################
# gives all combinations of isoenzyme indices in an ecFBA model
# (in a bit complicated way)
#
# input parameters:
#   - fba : FBA model struct with all important information
#
# output parameters:
#   - indices : array with index arrays of enzymes that catalyse the same 
#               reaction
#
################################################################################
function getIsoenzymes(fba::FBAModel)::Array{Array{Int, 1}, 1}

    # extract important content of model
    reactionNames = fba.reactionNames
    componentNames = fba.componentNames
    enzymeNames = fba.enzymeNames
    fluxes = fba.model[:fluxes]
    stochiometry = fba.model[:stochiometry]

    # find the pseudo arm components (since those are connected to isoenzymes),
    # exclude the reversible ones, since they will create doubles
    armComponents = findall(x -> contains(x, "pmet_") && !contains(x, "REV"), 
                            componentNames)
    
    # go through the pseudo arm components and find isoenzymes  
    enzymeCombinations = Array{Array{Int, 1}, 1}(undef, length(armComponents))
    for c in armComponents
    
        # get all the fluxes that the arm component is involved in
        involvedFluxIdx = findall(x -> x != 0, 
                                  normalized_coefficient.(stochiometry[c], fluxes))
        involvedReactions = reactionNames[involvedFluxIdx]

        # pick only the ones that are not the arm but the isoenzyme reactions
        isoenzymeReactions = involvedFluxIdx[.!contains.(involvedReactions, "(arm")]
    
        # go through the found reactions, find the enzymes used in the reaction and
        # map them to the indices of enzyme variables in the optimisation problem
        # (only the arm component and the isoenzymes have negative stochiometry
        # in those reactions)
        combination = Array{Int, 1}(undef, 0)
        enzymeComponents = contains.(componentNames, "prot_") .& 
                           .!contains.(componentNames, "_pool")
        for i in isoenzymeReactions
            tmpEnzymes = (normalized_coefficient.(stochiometry, fluxes[i]) .< 0) .&
                         enzymeComponents
            tmpEnzymeNames = map(x -> x[6:end], 
                                componentNames[findall(x -> x == true, tmpEnzymes)])
            tmpEnzymeIdx = map(y -> findfirst(x -> x == y, enzymeNames[:, 1]), 
                                tmpEnzymeNames)
            combination = [combination; tmpEnzymeIdx]
        end
    
        # in case of enzyme complexes 
        nIsoenzymes = length(combination)
        relevant = trues(nIsoenzymes)

        # take only the ones that differ between the reactions, 
        # since those are the once that replace each other
        for j = 1 : nIsoenzymes
            if sum(combination[j] .== combination) > 1
                relevant[j] = false
            end
        end
    
        # if only one enzyme remains, or the combination aready exists, 
        # don't add it
        idx = findfirst(x -> x == c, armComponents)
        combination = combination[relevant]
        if (length(combination) > 1) .& (isempty(findall(x -> x == combination, 
                                         enzymeCombinations[1:idx-1])))
            enzymeCombinations[idx] = combination
        else
            enzymeCombinations[idx] = []
        end
    end

    # remove empty elements
    enzymeCombinations = enzymeCombinations[.!isempty.(enzymeCombinations)]

    # return the combination of isoenzymes that belong to the same reaction
    return enzymeCombinations
end


################################################################################
# maps an objective (as string) to an array of coefficients c such that 
# c^T * v becomes the objective value
#
# OBS: this is specific to the model we use
#
# input parameters:
#   - string : describes the objective, e.g. "growth", "glucoseUptake", 
#              "ATPProduction", "NADHProduction", "allFluxes", "ATPHydrolysis"
#
# output parameters:
#   - coefficients : array with coefficients for all fluxes in the objective 
#                    function
#
################################################################################
function mapObjectiveToCoefficients(objective::String)::Array{Float64, 1}

    coefficients = zeros(Float64, 234)

    if objective == "growth"
        coefficients[40] = 1.0
    elseif objective == "glucoseUptake"
        coefficients[7] = 1.0
    elseif objective == "ATPProduction"
        idx = [72, 167, 80, 187, 126, 140, 154]
        coefficients[idx] = [1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0]    
    elseif objective == "NADHProduction"
        idx = [68, 163, 113, 184, 109, 122, 137, 139, 142, 192]
        coefficients[idx] = [1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0,-1.0]
    elseif objective == "allFluxes"
        coefficients = ones(Float64, 234)
    elseif objective == "ATPHydrolysis"
        coefficients[16] = 1.0
    elseif objective == "NGAM"
        coefficients[234] = 1.0
    end

    return coefficients
end