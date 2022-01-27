################################################################################
################################################################################
# InitialiseODEModel.jl includes functions for initialising the damage
# accumulation ODE model:
#   - odeModel = initialiseODEModel(growth, damageFormation, damageRepair, P, D, 
#                                   mass)
################################################################################
################################################################################


################################################################################
# define the ODE struct
################################################################################
mutable struct ODEModel
    fixedParameters::Array{Float64, 1}
    variableParameters::Array{Float64, 1}
    state::Array{Float64, 1}
end


################################################################################
# ODE model initialisation 
#
# input parameters:
#   - growth : growth rate
#   - damageFormation : damage formation rate
#   - damageRepair : damage repair rate
#   - intacts : fraction of dry mass coming from intact proteins
#   - damaged : fraction of dry mass coming from damaged proteins
#   - mass : dry mass
#
# output parameters:
#   - odeModel : ODE model struct with all important information
#
################################################################################
function initialiseODEModel(growth::Float64, damageFormation::Float64, 
                            damageRepair::Float64, intacts::Float64, 
                            damaged::Float64, mass::Float64)::ODEModel

    fixedParameters = [damageFormation, damageRepair, 0.0]
    variableParameters = [0.0, 0.0, growth]
    state = [intacts, damaged, mass]

    return ODEModel(fixedParameters, variableParameters, state)

end