################################################################################
################################################################################
# Cell.jl includes functions for simulating a cell's lifespan with the 
# integrated model (ecFBA for the metabolism, Boolean for the signalling, ODE
# for damage accumulation):
#   - cell = createCell(fbaPath, booleanSpeciesPath, booleanRulesPath, 
#                       TFpath, maxGrowthRate, damageFormation, damageRepair, 
#                       P, D, mass, signallingThresholds)
#   - relevantRefs = getRelevantReferences(fba, boolean)
#   - outputs = simulateLife!(cell, initialMass, timestep, maximalTimesteps, 
#                             regulationFactor, firstObjective, firstFlexibility, 
#                             secondObjective, secondFlexibility)
################################################################################
################################################################################


################################################################################
# cell struct
################################################################################
mutable struct Cell 
    fba::FBAModel
    boolean::BooleanModel
    tfTargets::Array{typejoin(Expr, Array{Int, 1}), 2}
    latestBooleanInputs::Array{Array{Bool, 1}, 1}
    ode::ODEModel
    signallingThresholds::Array{Float64, 1}
    relevantRefs::Dict{String, typejoin(VariableRef, ConstraintRef, BooleanComponent)}
end


################################################################################
# initialise a Cell struct including the three models as well as some parameters
# and references to relevant model parts
#
# input parameters:
#   - fbaPath : path to the .mat file of the ecFBA model
#   - booleanSpeciesPath : path to the .txt file specifying the components in the 
#                          Boolean model
#   - booleanRulesPath : path to the .txt file specifying the rules in the 
#                        Boolean model
#   - TFpath : path to the .txt file specifying the gene regulation of enzymes 
#              in the ecFBA model triggered by transcription factors in the 
#              Boolean model
#   - nDelay : indicates how many time steps it takes until the a signal actually
#              affects the metabolic network
#   - maxGrowthRate : maximal growth rate allowed in the ecFBA model
#   - damageFormation : rate at which damage forms, excluding the damage caused 
#                       by the metabolic network
#   - damageRepair : rate at which damage is repaired
#   - P : intact protein fraction
#   - D : damaged protein fraction
#   - mass : cell mass
#   - signallingThresholds : threshold values for the switching on of glucose, 
#                            hydrogen peroxide and Trx in the Boolean model
#   - firstObjective : Objective struct to describe the first objective
#   - secondObjective : Objective struct to describe the second objective
#                      
#
# output parameters:
#   - cell : Cell struct with all relevant information
#
################################################################################
function createCell(fbaPath::String, booleanSpeciesPath::String, 
                    booleanRulesPath::String, TFpath::String, nDelay::Int,
                    maxGrowthRate::Float64, damageFormation::Float64, 
                    damageRepair::Float64, P::Float64, D::Float64, mass::Float64,
                    signallingThresholds::Array{Float64, 1}, firstObjective::Objective, 
                    secondObjective::Objective)::Cell

    # initialise ODE model
    ode = initialiseODEModel(maxGrowthRate, damageFormation, damageRepair,
                             P, D, mass)

    # initialise ecFBA model
    fba = initialiseJuMPModel(fbaPath)
    prepareLifespanExperiment!(fba, maxGrowthRate)

    # add four coonstraints that correspond to bounding the objective function
    c = firstObjective.coefficients
    @constraint(fba.model, objective1Lb, c' * fba.model[:fluxes] >= 0)
    @constraint(fba.model, objective1Ub, c' * fba.model[:fluxes] <= 1000)

    c = secondObjective.coefficients
    @constraint(fba.model, objective2Lb, c' * fba.model[:fluxes] >= 0)
    @constraint(fba.model, objective2Ub, c' * fba.model[:fluxes] <= 1000)

    # initialise Boolean model
    boolean = initialiseBooleanModel(booleanSpeciesPath, booleanRulesPath)
    targets = readTranscriptionalTargets(TFpath, fba.enzymeNames)

    # extract and save relevant reactions and components that are used very often
    relevantRefs = getRelevantReferences(fba, boolean)

    # use the current Boolean model for the regulation in first few timesteps
    booleanInputs = [relevantRefs["exGlc"].present, 
                     relevantRefs["H2O2"].present,
                     relevantRefs["Trx1_2"].present]
    latestBooleanInputs = [booleanInputs for i = 1 : nDelay + 1]

    # update the GAM component in the FBA according to the maximal growth rate
    updateGAMValue!(relevantRefs["Growth"], 
                    relevantRefs["Maintainance for growth [cytosol]"], 
                    maxGrowthRate)

    # initialise and return cell
    return Cell(fba, boolean, targets, latestBooleanInputs, ode, 
                signallingThresholds, relevantRefs)
end


################################################################################
# map relevant reaction names in an ecFBA and Boolean model to the variables and 
# contraints in the model, in order to facilitate their usage
#
# input parameters:
#   - fba : FBA model struct with all important information
#   - boolean : Boolean model struct with all relevant information
#
# output parameters:
#   - relevantRefs : dictionary that translates each reaction or component name
#                    to the actual variable in the models
# 
################################################################################
function getRelevantReferences(fba::FBAModel, boolean::BooleanModel)::Dict{String, 
                               typejoin(VariableRef, ConstraintRef, BooleanComponent)}

    # extract relevant FBA reactions
    reactionsFBA = ["Growth", "NGAM", "Uptake of glucose", 
                    "Production of ethanol",
                    "Production of ethanol (reversible)",
                    "protein damage exchange (cytosol)",
                    "protein damage exchange (mitochondria)"]
    refsReactions = fba.model[:fluxes][map(x -> findfirst(y -> y == x, 
                              fba.reactionNames), reactionsFBA)]

    # extract relevant FBA enzymes
    enzymesFBA = ["draw_prot_P22217", "draw_prot_P22803"]
    refsEnzymes = fba.model[:enzymes][map(x -> findfirst(y -> y == x, fba.reactionNames), 
                            enzymesFBA) .- length(fba.model[:fluxes])]

    # extract relevant FBA components
    componentsFBA = ["Maintainance for growth [cytosol]"]
    refsComponents = fba.model[:stochiometry][map(x -> findfirst(y -> y == x, 
                               fba.componentNames), componentsFBA)]

    # extract relevant Boolean components for signalling
    componentsBool = ["exGlc", "H2O2", "Trx1_2"]
    refsBool = boolean.components[map(x -> findfirst(y -> y.name == x, 
                                  boolean.components), componentsBool)]

    # create dictionary with the important variable and constraint references 
    allRefs = [refsReactions; refsEnzymes; refsComponents; refsBool]
    allNames = [reactionsFBA; enzymesFBA; componentsFBA; componentsBool]
    relevantRefs = Dict(allNames .=> allRefs)

    return relevantRefs
end 


################################################################################
# simulate the life of a cell with the integrated model (ecFBA for the 
# metabolism, Boolean for the signalling, ODE for damage accumulation) until 
# the maximal number of timesteps are reached or the cell died,
# with one or two successive objectives (it is always the parsimonious solution
# that is chosen)
# 
# a cells life is divided into three phases:
#   1. maximal growth phase (via fermentation)
#   2. switch to respiration due to increased ATP demand for repair
#   3. ethanol respiration 
#
#
# input parameters:
#   - cell : Cell struct with all relevant information
#   - minimalMass : minimal mass of the mother cell (important for division)
#   - sizeProportion : fraction of biomass that remains in the mother cell at 
#                      cell divison
#   - retention : fraction of damage in the mother that is caused by 
#                 damage retention
#   - timestep : time of the iteration
#   - maximalTimesteps : maximal number of timesteps before the simulation stops
#   - regulationFactor : factor determining how much the Boolean signalling layer
#                        affects the enzyme usages in the ecFBA model
#   - firstObjective : Objective struct to describe the first objective
#   - secondObjective : Objective struct to describe the second objective
#
# output parameters:
#   - status : indicates if the cell is still alive, or the reason of death
#   - rls : replicative lifespan of the cell
#   - averageGenerationTime : average time between divisions
#   - phaseSplit : times in, number of divisions during and damage levels at the 
#                  end of the three phases described above
# 
################################################################################
function simulateLife!(cell::Cell, minimalMass::Float64, sizeProportion::Float64, 
                       retention::Float64, timestep::Float64, 
                       maximalTimesteps::Int, regulationFactor::Float64,
                       firstObjective::Objective, secondObjective::Objective
                       )::Tuple{String, Int64, Float64, Array{Array{Float64, 1}, 1}}

    
    ############################################################################
    # INITIALISATION
    ############################################################################
    # initialise all counts
    status = "alive"
    time = 0.0
    rls = 0
    divisionTimes = [0.0]
    it = 1
    phaseSplit = [[0.0, 0.0, NaN] for i = 1 : 3]
    currentPhase = "phase 1"
    phaseShift = [false, false]

    # extract relevant model parts that are used over and over again
    fba = cell.fba
    model = fba.model
    boolean = cell.boolean
    ode = cell.ode

    # extract important references to the models
    refs = cell.relevantRefs
    glucoseUptakeFlux = refs["Uptake of glucose"]
    ethanolProductionFlux = refs["Production of ethanol"]
    ethanolUptakeFlux = refs["Production of ethanol (reversible)"]
    growthFlux = refs["Growth"]
    GAMComponent = refs["Maintainance for growth [cytosol]"]
    NGAM = refs["NGAM"]
    damageFluxes = [refs["protein damage exchange (cytosol)"],
                    refs["protein damage exchange (mitochondria)"]]
    trxFluxes = [refs["draw_prot_P22217"],
                 refs["draw_prot_P22803"]]
    booleanGlucose = refs["exGlc"]
    booleanPeroxoide = refs["H2O2"]
    booleanTrx = refs["Trx1_2"]

    # save initial enzyme, growth and glucose bounds
    fluxLowerBounds = lower_bound.(model[:fluxes])
    fluxUpperBounds = upper_bound.(model[:fluxes])
    enzymeLowerBounds = lower_bound.(model[:enzymes])
    enzymeUpperBounds = upper_bound.(model[:enzymes])

    # save initial bound on objective
    objectiveBounds = normalized_rhs.([model[:objective1Lb], model[:objective1Ub], 
                                       model[:objective2Lb], model[:objective2Ub]])

    # get maximal growth rate and NGAM 
    maxGrowthRate = upper_bound(growthFlux)
    initialGrowth = maxGrowthRate
    maxNGAM = upper_bound(NGAM)


    ############################################################################
    # SIMULATE TIMESTEP BY TIMESTEP
    ############################################################################
    for n = 1 : maximalTimesteps
    
        ########################################################################
        # UPDATE ecFBA PARAMETERS
    
        # limit the enzyme pool according to current ODE state
        P = ode.state[1]
        set_upper_bound(model[:pool], P * 0.1799 * 0.4592)

        # set NGAM linearly dependent on the current amount of damage in the cell
        # between 0 and a maximal value (set as upper bound, from Lu et al., 2019)
        # plus a baseline
        baseline = 0.0
        D = ode.state[2]
        set_lower_bound(NGAM, baseline + D * maxNGAM / (P + D))

        ########################################################################
        # SOLVE ecFBA MODEL FOR THE FIRST OBJECTIVE
        FBA!(model, firstObjective, 1) 

        # stop if the FBA cannot solve at all anymore (no possible growth rate)
        if termination_status(model) != MOI.OPTIMAL 
            status = "dead (regular)"
            break
        end

        # restrict the objective value with some flexibility
        objValue = firstObjective.coefficients' * value.(model[:fluxes])
        if firstObjective.sense == MOI.MAX_SENSE 
            set_normalized_rhs(model[:objective1Lb], objValue * 
                                                     (1 - firstObjective.flexibility))
        elseif firstObjective.sense == MOI.MIN_SENSE 
            set_normalized_rhs(model[:objective1Ub], objValue * 
                                                     (1 + firstObjective.flexibility))
        end

        # update the GAM to the current growth rate
        updateGAMValue!(growthFlux, GAMComponent, value(growthFlux))

        ########################################################################
        # SOLVE A SECOND TIME IN CASE OF A SECOND OBJECTIVE
        if secondObjective.description != "" 

            FBA!(model, secondObjective, 2) 

            objValue = secondObjective.coefficients' * value.(model[:fluxes])
            if secondObjective.sense == MOI.MAX_SENSE 
                set_normalized_rhs(model[:objective2Lb], objValue * 
                                                         (1 - secondObjective.flexibility))
            elseif secondObjective.sense == MOI.MIN_SENSE 
                set_normalized_rhs(model[:objective2Ub], objValue * 
                                                         (1 + secondObjective.flexibility))
            end
        end

        ########################################################################
        # FIND THE LOGICAL STEADY STATE OF THE BOOLEAN MODEL AND REGULATE THE
        # ecFBA ACCORDINGLY

        # update the Boolean model input with the ones from nDelay timesteps
        # ago, i.e. the first element in the latestBooleanInputs array
        tmpBooleanInputs = cell.latestBooleanInputs[1]
        refs["exGlc"].present = tmpBooleanInputs[1]
        refs["H2O2"].present = tmpBooleanInputs[2]
        refs["Trx1_2"].present = tmpBooleanInputs[3]

        # run Boolean model 
        runBooleanModel!(boolean)

        # calculate ranks for gene expression given the transcription factor 
        # activity
        ranks = getExpressionRank(boolean.components, fba.enzymeNames, 
                                  cell.tfTargets)

        # regulate FBA according to the ranks from the Boolean model
        regulateFBA!(model, ranks, regulationFactor)

        ########################################################################
        # SOLVE THE REGULATED ecFBA MODEL
        if secondObjective.description == "" 
            FBA!(model, firstObjective, 1) 
        else
            FBA!(model, secondObjective, 2) 
        end

        # stop if the FBA cannot solve anymore
        if termination_status(model) != MOI.OPTIMAL 
            status = "dead (regulation)"
            break
        end

        ########################################################################
        # SOLVE ODE MODEL WITH PARAMETERS TAKEN FROM THE OPTIMAL ecFBA SOLUTION

        # extract the growth and damage fluxes of the optimal solution
        damage = sum(value.(damageFluxes))
        repair = 0.0
        growth = value(growthFlux)

        # set variable rate parameters for the ODE model accordingly
        ode.variableParameters = [damage, repair, growth]

        # update ODE for one timestep
        solveODE!(ode, timestep)
    
        ########################################################################
        # CELL DIVISION
        # if the cell has build up enough biomass (i.e. so that a new daughter
        # has at least 1/s - 1 of the mother cell mass)
        # and the cell has a large enough fration of intact proteins (so that 
        # retention does not cause a negative fraction of intact proteins)
        # (respect precision of solver)
        P, D, mass = ode.state
        bigEnough = mass - 1 / sizeProportion * minimalMass >= -precision
        enoughIntacts = P - D * retention >= -precision
    
        # the cell divides 
        if bigEnough && enoughIntacts
            division!(ode, sizeProportion, retention)
            rls += 1
            push!(divisionTimes, time)
        end

        ########################################################################
        # UPDATE THE BOOLEAN MODEL ACCORDING TO THE OPTIMAL ecFBA SOLUTION
        # update cell's signalling according to solution
        triggerSignalling!([glucoseUptakeFlux], cell.signallingThresholds[1], 
                booleanGlucose)
        triggerSignalling!(damageFluxes, cell.signallingThresholds[2], 
                        booleanPeroxoide)
        triggerSignalling!(trxFluxes, cell.signallingThresholds[3], booleanTrx)

        # update current Boolean model to be used in nDelay time steps,
        # i.e. the last element in the latestBooleanInputs array
        cell.latestBooleanInputs[end] = [refs["exGlc"].present, 
                                         refs["H2O2"].present,
                                         refs["Trx1_2"].present]

        # and shift all elements by one for the next time step
        cell.latestBooleanInputs[1:end-1] = cell.latestBooleanInputs[2:end]

        ########################################################################
        # DECIDE WHICH PHASE THE CELL IS IN
        # distinguish between phases (respect precision of solver)
        it == 1 ? initialGrowth = growth : nothing
        highestGrowth = growth - initialGrowth * 0.95 >= -precision
        ethanolUptake = (value(ethanolUptakeFlux) - value(ethanolProductionFlux)) > 
                         precision

        # phase 1
        # update time, rls and damage if the growth is still above the maximal 
        # values with some flexibility (and phase 2 has not been reached before)
        if highestGrowth && !ethanolUptake && phaseShift[1] == false
            currentPhase = "phase 1"
            phaseSplit[1][1] += timestep
            phaseSplit[1][2] = rls
            phaseSplit[1][3] = cell.ode.state[2]
        # phase 2
        # update time, rls and damage if the cell is switching to respiration
        # (and phase 3 has not been reached before)
        elseif !highestGrowth && !ethanolUptake && phaseShift[2] == false
            currentPhase = "phase 2"
            phaseShift[1] = true
            phaseSplit[2][1] += timestep
            phaseSplit[2][2] = rls - phaseSplit[1][2] 
            phaseSplit[2][3] = cell.ode.state[2]
        # phase 3
        # update time, rls and damage if there is ethanol respiration
        elseif ethanolUptake
            currentPhase = "phase 3"
            phaseShift[2] = true
            phaseSplit[3][1] += timestep
            phaseSplit[3][2] = rls - phaseSplit[1][2] - phaseSplit[2][2] 
            phaseSplit[3][3] = cell.ode.state[2]
        end

        ########################################################################
        # UPDATE ecFBA PARAMETERS FOR NEXT TIMESTEP
        # unconstrain enzymes and fluxes again
        set_lower_bound.(model[:fluxes], fluxLowerBounds)
        set_upper_bound.(model[:fluxes], fluxUpperBounds)
        set_lower_bound.(model[:enzymes], enzymeLowerBounds)
        set_upper_bound.(model[:enzymes], enzymeUpperBounds)

        # unconstrain bounds on objectives again
        set_normalized_rhs.([model[:objective1Lb], model[:objective1Ub], 
                             model[:objective2Lb], model[:objective2Ub]], 
                            objectiveBounds)

        # update time and iteration
        time += timestep
        it += 1
        
    end 
    
    ########################################################################
    # GET AVERAGE GENERATION TIME AND RETURN RESULTS
    if rls > 0
        averageGenerationTime = mean(divisionTimes[2:end-1] .- 
                                     divisionTimes[1:end-2])
    else
        averageGenerationTime = 0.0
    end

    ########################################################################
    # RETURN
    output = (status, rls, averageGenerationTime, phaseSplit)
    return output

end