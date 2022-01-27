################################################################################
################################################################################
# InitialiseBooleanModel.jl includes functions for initialising a Boolean model:
#   - outputModel = initialiseBooleanModel(filepathSpecies, filepathRules)
################################################################################
################################################################################


################################################################################
# define the model structs
################################################################################
mutable struct BooleanComponent
    id::Int
    name::String
    systematicName::String
    pathway::String
    type::String
    present::typejoin(Nothing, Bool)
    phosphorylated::typejoin(Nothing, Bool)
    oxidised::typejoin(Nothing, Bool)
    active::typejoin(Nothing, Bool)
end

struct BooleanRule
    id::Int
    type::String
    description::String
    condition::Expr
    update::Expr
    alternative::typejoin(Expr, Symbol)
    reference::String
    active::Bool
end

struct BooleanModel
    components::Array{BooleanComponent, 1}
    nComponents::Int
    rules::Array{BooleanRule, 1}
    nRules::Int
end


################################################################################
# Boolean model initialisation from two model files
#
# input parameters:
#   - filepathSpecies : string with the path to the txt file with all species
#   - filepathRules : string with the path to the txt file with all rules
#
# output parameters:
#   - outputModel : Boolean model struct containing all rules and components with
#                   their initial state
#
################################################################################
function initialiseBooleanModel(filepathSpecies::String, 
                                filepathRules::String)::BooleanModel

    # extract species data
    data = readdlm(filepathSpecies,'\t', String, '\n', skipstart = 1)
    nComponents = size(data, 1)

    # parse data to right format
    ids = tryparse.(Int64, data[:, 1])
    names = data[:, 2]
    presence = tryparse.(Bool, data[:, 3])
    phosphorylation = tryparse.(Bool, data[:, 4])
    oxidation = tryparse.(Bool, data[:, 5])
    activity = tryparse.(Bool, data[:, 6])
    pathways = data[:, 7]
    types = data[:, 8]
    systematicNames = data[:, 9]

    # initialise vector with components for Boolean model and a dictionary 
    # that maps the name to the component index
    components = Array{BooleanComponent, 1}(undef, nComponents)

    # fill with data from file
    for i = 1 : nComponents
        components[ids[i]] = BooleanComponent(ids[i], names[i], systematicNames[i], 
                                              pathways[i], types[i], presence[i], 
                                              phosphorylation[i], oxidation[i], 
                                              activity[i])
    end

    # extract rules
    data = readdlm(filepathRules,'\t', String, '\n', skipstart = 1)
    nRules = size(data, 1)

    # parse data to right format
    ids = tryparse.(Int64, data[:, 1])
    types = data[:, 2]
    descriptions = data[:, 3]
    ifStatements = Meta.parse.(data[:, 4])
    thenStatements = Meta.parse.(data[:, 5])
    elseStatements = Meta.parse.(data[:, 6])
    references = data[:, 7]
    active = tryparse.(Bool, data[:, 8])

    # initialise vector with rules for Boolean model
    rules = Array{BooleanRule, 1}(undef, nRules)

    # fill with data from file
    for i = 1:nRules
        rules[ids[i]] = BooleanRule(ids[i], types[i], descriptions[i], ifStatements[i], 
                                    thenStatements[i], elseStatements[i], references[i], 
                                    active[i])
    end

    # return
    return BooleanModel(components, nComponents, rules, nRules)
end