################################################################################
################################################################################
# Boolean.jl includes functions for running a Boolean model:
#   - iterateBooleanModel!(boolean)
#   - knockout!(proteins, booleanComponents)
#   - iteration, steadyState = runBooleanModel!(boolean, maxIteration)
################################################################################
################################################################################


################################################################################
# run one synchronous iteration of the Boolean rules and updates the states
#
# input parameters:
#   - boolean : Boolean model struct with all relevant information
#
################################################################################
function iterateBooleanModel!(boolean::BooleanModel)::Nothing

    # extract the model's components to work with
    tmpComponents = boolean.components

    # create the right namespace for the conditions with each component in the 
    # model having actually the name of the protein as a variable in the code 
    # here, such that the rules can be directly executed
    names = map(x -> x.name, tmpComponents)
    [@eval $(Symbol(names[i])) = $tmpComponents[$i] for i in 1:boolean.nComponents]
    
    # create boolean array of which rules apply and which not
    conditions = [eval(boolean.rules[i].condition) for i in 1:boolean.nRules]

    # update all components according to the set rules in the model, 
    # but only if the respective rule is active
    for i = 1 : boolean.nRules
        if boolean.rules[i].active 
            conditions[i] ? eval(boolean.rules[i].update) : 
                            eval(boolean.rules[i].alternative)
        end
    end

    # return
    return nothing
end


################################################################################
# knockout proteins in a Boolean model
#
# input parameters:
#   - proteins : string specifiying the proteins that should be knocked out, 
#                separated by a "+" in the string
#   - booleanComponents : components of a Boolean model struct
#
################################################################################
function knockout!(proteins::String, 
                   booleanComponents::Array{BooleanComponent, 1})::Nothing

    # get the names of the proteins to find the right index for the knockouts
    names = map(x -> x.name, booleanComponents)

    if proteins != ""
        # divide the string with the proteins to get the individual ones
        individualProteins = split(proteins, "+")

        for protein in individualProteins
            # find the protein index
            proteinIdx = findfirst(x -> x == protein, names)

            # set all properties to false
            booleanComponents[proteinIdx].present = false
            booleanComponents[proteinIdx].phosphorylated = false
            booleanComponents[proteinIdx].oxidised = false
            booleanComponents[proteinIdx].active = false
        end
    end

    # return 
    return nothing
end


################################################################################
# run Boolean model until a steady state or a maximal number of iterations
# is reached
#
# input parameters:
#   - boolean : Boolean model struct with all relevant information
#   - maxIteration : maximal number of iterations
#
# output parameters:
#   - steadyState : indicates if a steady state is reached
#   - iteration : the number of iterations needed
#
################################################################################
function runBooleanModel!(boolean::BooleanModel, 
                          maxIteration::Int = 100)::Tuple{Bool, Int}

    # intitalise the steady state
    steadyState = false

    # save the current state 
    components = boolean.components
    state = [map(x -> x.present, components), 
             map(x -> x.phosphorylated, components),
             map(x -> x.oxidised, components),
             map(x -> x.active, components)]

    # iterate until steady state is reached
    iteration = 0
    while !steadyState && iteration < maxIteration
        # iterate one step
        iteration += 1
        iterateBooleanModel!(boolean)
    
        # evaluate if steady state is reached
        updatedState = [map(x -> x.present, components), 
                        map(x -> x.phosphorylated, components),
                        map(x -> x.oxidised, components),
                        map(x -> x.active, components)]
        steadyState = isequal(state, updatedState)

        # update the state
        state = updatedState
    end

    # return 
    return steadyState, iteration
end