################################################################################
################################################################################
# Integration.jl includes functions for connecting the FBA with the Boolean and
# the ODE model:
#   - triggerSignalling!(fbaFluxes, threshold, booleanComponent)
#   - targets = readTranscriptionalTargets(TFfilepath, fbaEnzymes)
#   - ranks = getExpressionRank(booleanComponents, fbaEnzymes, targets)
#   - regulateFBA!(model, ranks, regulationFactor)
################################################################################
################################################################################


################################################################################
# interface between the FBA and the Boolean part to trigger signalling
#
# input parameters:
#   - fluxes : fluxes in the FBA model that are used to switch on or off
#                 the signalling
#   - threshold : determines from which value of the sum of the fluxes the 
#                 feedback happens
#   - booleanComponent : component in the Boolean model that is switched on or 
#                        off in response to the fluxes
#
################################################################################
function triggerSignalling!(fluxes::Array{VariableRef, 1}, 
                            threshold::Float64, 
                            booleanComponent::BooleanComponent)::Nothing

    # if the given flux (or the sum of fluxes) is above a certain threshold
    # set the corresponding Boolean component to present
    if sum(value.(fluxes)) > threshold
        booleanComponent.present = true
    else
        booleanComponent.present = false
    end

    # return
    return nothing
end


################################################################################
# read in targets of transcription factors
#
# input parameters:
#   - TFfilepath : path to the collection of transcription factor targets
#   - fbaEnzymes : enzyme attribute of FBA model
#
# output parameters:
#   - targets : 2d array with conditions in the Boolean model and the resulting 
#               indices of enzymes in the FBA model that are up or downregulated
#
################################################################################
function readTranscriptionalTargets(TFfilepath::String, 
                                    fbaEnzymes::Array{String, 2})::Array{Any, 2}

    # read file
    input = readdlm(TFfilepath,'\t', String, '\n', skipstart = 1)
    nTF = size(input, 1)

    # extract transcription factors (can be several for coregulation)
    tfs = split.(input[:, 1], " && ")

    # and initialise up and downregulated genes
    conditions = Array{Expr, 1}(undef, nTF)
    positives = Array{Array{Int, 1}, 1}(undef, nTF)
    negatives = Array{Array{Int, 1}, 1}(undef, nTF)

    # go through them and map them to the enzyme indices in the FBA model
    for i = 1 : nTF

        # get the condition for the specific transcription factor to be 
        # active
        conditions[i] = Meta.parse.(join(tfs[i] .* ".active", " && "))

        # up regulation
        upRegulated = split(input[i, 2], ",")
        if !isempty(upRegulated)
            idx = map(x -> findfirst(y -> y == x, fbaEnzymes[:, 2]), 
                      upRegulated)
            positives[i] = idx[.!isnothing.(idx)]
        else
            positives[i] = []
        end

        # down regulation
        downRegulated = split(input[i, 3], ",")
        if !isempty(downRegulated)
            idx = map(x -> findfirst(y -> y == x, fbaEnzymes[:, 2]), 
                      downRegulated)
            negatives[i] = idx[.!isnothing.(idx)]
        else
            negatives[i] = []
        end

    end
    
    # return
    return [conditions positives negatives]
end


################################################################################
# calculate the rank of each gene expression based on the Boolean model and 
# the target genes for the included transcription factors, as well as 
# posttranscriptional modifications 
# (from Österberg et al., PLOS Computational Biology 2021)
#
# input parameters:
#   - booleanComponents : components of a Boolean model struct
#   - fbaEnzymes : enzyme attribute of FBA model
#   - targets : 2d array with conditions in the Boolean model and the resulting 
#               indices of enzymes in the FBA model that are up or downregulated
#
# output parameters:
#   - ranks : array with a rank for each enzyme in the FBA model
#
################################################################################
function getExpressionRank(booleanComponents::Array{BooleanComponent, 1}, 
                           fbaEnzymes::Array{String, 2}, 
                           targets::Array{Any, 2})::Array{Int, 1}

    # extract necessary information from FBA model
    nEnzymes = size(fbaEnzymes, 1)
    ranks = zeros(Int, nEnzymes)

    # get number of conditions to be checked
    nConditions = size(targets, 1)

    # create correct namespace in Boolean model
    nComponents = length(booleanComponents)
    names = map(x -> x.name, booleanComponents)
    [@eval $(Symbol(names[i])) = $booleanComponents[$i] for i in 1:nComponents]

    # go through conditions in Boolean model and calculate ranks for each 
    # enzyme in the FBA model
    for i = 1 : nConditions
        
        # check condition stored in targets
        if eval(targets[i, 1])
            # add one for each enzyme that is upregulated
            ranks[targets[i, 2]] .+= 1
            # sustract one for each enzyme that is downregulated
            ranks[targets[i, 3]] .-= 1
        end

    end
    
    # return
    return ranks
end


################################################################################
# adapts the lower and upper bounds of enzyme usages in the FBA model according
# to gene ranks that are based on the activity of the respective transcription 
# factors in the Boolean model
# (from Österberg et al., PLOS Computational Biology 2021)
#
# additional information:
# when creating/modifying the FBA model it was made sure that:
#   - enzyme names in fba.enzymeNames and in fba.reactionsNames appear in 
#     the same order
#   - fba.reactionNames are corresponding to the fba.model variables fluxes,
#     enzymes and pool, in that order
#   - consequently the ranks (made with fba.enzymeNames) can be directly 
#     mapped to the enzyme variables in fba.model
#
# input parameters:
#   - model : JuMP optimisation model of the metabolic network
#   - ranks : array with a rank for each enzyme in the FBA model 
#             (from getExpressionRank())
#   - regulationFactor : determines the fraction of the enzyme variability that
#                        is added or reduced on the respective bounds
#
################################################################################
function regulateFBA!(model::Model, ranks::Array{Int, 1}, 
                      regulationFactor::Float64)::Nothing

    # find which enzymes are regulated (rank != 0)
    regulatedIdx = findall(x -> x != 0, ranks)

    if !isempty(regulatedIdx)

        # distinguish between up and down regulation
        upregulated = findall(x -> x > 0, ranks[regulatedIdx])
        downregulated = findall(x -> x < 0, ranks[regulatedIdx])

        # and extract the enzyme variable 
        enzymes = model[:enzymes]
        optimalUsages = value.(enzymes)

        # get minimal and maximal usage of those enzymes 
        # for the given (constraint) FBA 
        minUsage, maxUsage = calculateEVA(model, enzymes[regulatedIdx])

        # use optimal solution in case EVA was not optimal
        # (sometimes it can happen that there are numerical errors ...)
        noResult = isnan.(minUsage)
        minUsage[noResult] = optimalUsages[regulatedIdx[noResult]]
        noResult = isnan.(maxUsage)
        maxUsage[noResult] = optimalUsages[regulatedIdx[noResult]]

        # get the difference 
        delta = maxUsage .- minUsage

        # negative rank means downregulation of the upper bound of that enzyme by
        # a factor 
        set_upper_bound.(enzymes[regulatedIdx[downregulated]], 
                         maxUsage[downregulated] .- 
                         (delta[downregulated] .* regulationFactor))

        # positive rank means upregulation of the lower bound of that enzyme by
        # a factor 
        set_lower_bound.(enzymes[regulatedIdx[upregulated]], 
                         minUsage[upregulated] .+ 
                         (delta[upregulated] .* regulationFactor))
    end
    
    # return
    return nothing
end