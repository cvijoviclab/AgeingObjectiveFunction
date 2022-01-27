################################################################################
################################################################################
# FBA.jl includes useful functions for performing a (parsimonious) flux balance
# analysis:
#   - updateGAMValue!(growthReaction, GAMComponent, growthRate)
#   - FBA!(model, objective, parsimonious, stage)
#   - minUsage, maxUsage = calculateEVA(model, enzymes)
#   - enzymesToComponents, componentsToFluxes = mapEnzymeIdx(fba)
################################################################################
################################################################################


# define objective struct for the optimisations
struct Objective
    description::String
    coefficients::Array{Float64, 1}
    sense:: MOI.OptimizationSense
    flexibility::Float64
    parsimonious::Bool
end


################################################################################
# calculate and update GAM depending on the growth rate 
# according to Österberg et al., PLOS Computational Biology 2021 
#
# input parameters:
#   - growthReaction : reference to the variable in the JuMPmodel that 
#                      corresponds to the growth
#   - GAMComponent : reference to the mass balance constraint in the JuMPmodel  
#                    that corresponds to the GAM
#   - growthRate : current growth rate (obs: formula only valid for growth rates
#                  between 0.0 and 0.4)
#
################################################################################
function updateGAMValue!(growthReaction::VariableRef, GAMComponent::ConstraintRef, 
                         growthRate::Float64)::Nothing

    # linear change from 18 to 30 mmol / gDW h
    # for growth rates between 0.0 and 0.285
    if 0.0 <= growthRate <= 0.285
        GAM = 18 + growthRate * (30 - 18) / 0.285
        
        # update the reaction stochiometry accordingly
        set_normalized_coefficient(GAMComponent, growthReaction, - GAM)

    # linear change from 30 to 25 mmol / gDW h
    # for growth rates between 0.285 and 0.4
    elseif 0.285 < growthRate <= 0.4
        GAM = 30 + (growthRate - 0.285) * (25 - 30) / (0.4 - 0.285)

        # update the reaction stochiometry accordingly
        set_normalized_coefficient(GAMComponent, growthReaction, - GAM)
    end

    # return
    return nothing
end


################################################################################
# solve FBA for a JuMP model:
# either on optimisation or parsimonious solution, i.e.
# two successive optimisations
#   (1) optimise a given fluxes from input
#   (2) given the optimal solution of (1), minimize the enzyme usage and the
#       sum of all fluxes 
#
# input parameters:
#   - model : JuMP optimisation model of the metabolic network
#   - objective : Objective struct to describe the objective and the properties
#   - stage : indicates if it is the first (1) or the second (2) stage of the 
#             optimisation  
#
################################################################################
function FBA!(model::Model, objective::Objective, stage::Int = 1)::Nothing

    # perform first optimisation of parsimonious FBA
    @objective(model, objective.sense, objective.coefficients' * model[:fluxes])
    optimize!(model)

    # only if the first optimisation was successful, do the second
    if termination_status(model) == MOI.OPTIMAL && objective.parsimonious == true

        # get the right constraints to modify
        if stage == 1
            lbConstraint = model[:objective1Lb]
            ubConstraint = model[:objective1Ub]
        else
            lbConstraint = model[:objective2Lb]
            ubConstraint = model[:objective2Ub]
        end

        # save bounds for setting model back to input state
        tmpLowerBound = normalized_rhs(lbConstraint)
        tmpUpperBound = normalized_rhs(ubConstraint)

        # get the optimal value to set as a constraint in the second
        # optimisation
        optimalValue = objective_value(model)

        # keep optimal glucose uptake for the second optimisation
        # we assume that in an irreversible model all fluxes are positive
        set_normalized_rhs(lbConstraint, optimalValue * (1 - precision))
        set_normalized_rhs(ubConstraint, optimalValue * (1 + precision))

        # set objective value for second optimisation
        parsimoniousObjective = model[:pool] + sum(model[:fluxes])
        @objective(model, Min, parsimoniousObjective)

        # perform second optimisation of parsimonious FBA
        optimize!(model)

        #= error message if not feasible
        if termination_status(model) != MOI.OPTIMAL
        @printf("ERROR: parsimonious FBA with status %s.\n",check 
                termination_status(model))
        end =#

        # reset bounds
        set_normalized_rhs(lbConstraint, tmpLowerBound)
        set_normalized_rhs(ubConstraint, tmpUpperBound)

        #= error message if not feasible
        else
        @printf("ERROR: FBA with status %s.\n", termination_status(model)) =#
    end

    # return
    return nothing
end


################################################################################
# perform an enzyme variablility analysis (EVA) for the FBA model by minimising
# and maximising given enzyme usages of interest to see how sensitive a 
# previously constraint system is to changes the respective usages
# (for further details see Mahadevan et al., Metabolic Engeneering 2003, 
# Österberg et al., PLOS Computational Biology 2021)
#
# input parameters:
#   - model : JuMP optimisation model of the metabolic network
#   - enzymes : references to the enzymes that should be investigated with EVA
#
# output parameters:
#   - minUsage : minimal usages of given enzymes under the given conditions
#   - maxUsage : maximal usage of given enzymes under the given conditions
#
################################################################################
function calculateEVA(model::Model, enzymes::Array{VariableRef, 1} 
                      = Array{VariableRef, 1}(undef, 0))::Tuple{Array{Float64, 1}, 
                      Array{Float64, 1}}
    
    # if no indices are given check all fluxes
    isempty(enzymes) ? enzymes = model[:enzymes] : nothing

    # initialise output
    nEnzymes = length(enzymes)
    minUsage = [NaN for i = 1 : nEnzymes]
    maxUsage= [NaN for i = 1 : nEnzymes]
    
    # go through all enzymes and minimise and maximise the respective usage given
    # the input model (that should potentially have constraints on other fluxes)
    for i = 1 : nEnzymes

        # solve for minimal usage
        @objective(model, Min, enzymes[i])
        optimize!(model)
        
        # save if optimisation was successful
        if termination_status(model) == MOI.OPTIMAL
            minUsage[i] = value(enzymes[i])
        end
        
        # solve for maximal usage
        @objective(model, Max, enzymes[i])
        optimize!(model)
        
        # save if optimisation was successful
        if termination_status(model) == MOI.OPTIMAL
            maxUsage[i] = value(enzymes[i])
        end
    end
    
    # return
    return minUsage, maxUsage
end


################################################################################
# maps each enzyme to the components and reactions that it is involved 
#
# input parameters:
#   - fba : FBA model struct with all important information
#
# output parameters:
#   - enzymesToComponents : dictionary to map between the enzyme idx in the 
#                           enzymes variable/names and its idx in the components
#   - componentsToFluxes : dictionary to map between the enzyme idx in the 
#                          components to the fluxes/reactions that it is 
#                          involved in
#
################################################################################
function mapEnzymeIdx(fba::FBAModel)::Tuple{Dict{Int, Int}, Dict{Int, 
                                      Array{Int, 1}}}

    # extract relevant information from the FBA model
    model = fba.model
    fluxes = model[:fluxes]
    stochiometry = model[:stochiometry]
    components = fba.componentNames
    enzymes = fba.enzymeNames

    # initialise output dictionary
    enzymesToComponents = Dict{Int, Int}()
    componentsToFluxes = Dict{Int, Array{Int, 1}}()

    # go through enzymes and save the corresponding reaction indices 
    # (only in fluxes)
    for i = 1:size(enzymes, 1)
        nameInModel = enzymes[i, 1]
        idx = findfirst(x -> x == "prot_" * nameInModel, components)
        coefficients = normalized_coefficient.(stochiometry[idx], fluxes)
        fluxIndices = findall(x -> x != 0.0, coefficients)

        # add to dicitionary
        enzymesToComponents[i] = idx
        componentsToFluxes[idx] = fluxIndices
    end

    # return 
    return enzymesToComponents, componentsToFluxes
end
