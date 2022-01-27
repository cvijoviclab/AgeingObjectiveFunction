################################################################################
################################################################################
# DamageODE.jl includes functions for solving the damage accumulation ODE model:
#   - solveODE!(ode, time)
#   - daughter = division!(ode, sizeProportion, retention)
################################################################################
################################################################################


################################################################################
# solve ODE model
#
# dp/dt = -f * p + r * d
# dd/dt = f * p - r * d
#
# with eigenvalues (0, -(f+r))
# and eigenvectors (f/r, 1) and (-1, 1)
#
# dm/dt = g * m
# with simple exponential solution
#
# input parameters:
#   - ode : ODE model struct with all important information
#   - time : time until model should be solved
#
# output parameters:
#   - newState : ODE state at time defined as input
#
################################################################################
function solveODE!(ode::ODEModel, time::Float64)::Nothing

    f, r, g = ode.fixedParameters .+ ode.variableParameters
    p, d, m = ode.state
    
    # add f=r=0 as a extra case to avoid having 0 in the denominator
    # then dp/dt=dd/dt=0 and the solution does not change 
    # (is also limit of analytical solution, but the program complains)
    if f == 0 && r == 0
        pNew = p
        dNew = d
    # otherwise update the protein fractions according to analytical solution
    # (simplified eigenvalue/-vector solution)
    else
        pNew = 1 / (f + r) * (r * (p + d) - (d * r - p * f) * exp(-(f + r) * time))
        dNew = 1 / (f + r) * (f * (p + d) + (d * r - p * f) * exp(-(f + r) * time))
    end

    # calculate new mass with simple exponential growth
    mNew = m * exp(g * time)

    # update state
    ode.state = [pNew, dNew, mNew]

    # return
    return nothing
end


################################################################################
# define protein fractions and mass after division in the mother and the bud,
# according to asymmetry in size and additional retention mechanisms
#
# input parameters:
#   - ode : ODE model struct with all important information
#   - sizeProportion : fraction of mother cell compartment
#   - retention : retention factor
# 
# output parameters:
#   - mother : state of mother cell after division
#   - daughter : state of daughter cell after division
#
################################################################################
function division!(ode::ODEModel, sizeProportion::Float64, 
                   retention::Float64)::Array{Float64,1}

    p, d, m = ode.state

    # update the protein fractions according to retention
    retainedDamageFraction = retention * d
    motherP = p - retainedDamageFraction
    motherD = d + retainedDamageFraction
    daughterP = p + retainedDamageFraction
    daughterD = d - retainedDamageFraction

    # update the masses
    motherM = sizeProportion * m
    daughterM = (1 - sizeProportion) * m

    # define new mother and daughter states
    mother = [motherP, motherD, motherM]
    daughter = [daughterP, daughterD, daughterM]

    # update mother state and return daughter
    ode.state = mother
    return daughter
end