# FBA parameters
const fbaPath = "../../../ModelFiles/ecFBA/reducedEcYeast_modified_withDamage.mat"

# set the first FBA objective
const firstObjective = ["ATPProduction" for i = 1:5]
const firstSense = [MOI.MAX_SENSE for i = 1:5]
const firstFlexibility = collect(0.05:0.05:0.5)

# set the second FBA objective
const secondObjective = ["", "growth", "glucoseUptake", "NADHProduction", "NGAM"]
const secondSense = [MOI.MAX_SENSE, MOI.MAX_SENSE, MOI.MIN_SENSE, MOI.MIN_SENSE, MOI.MAX_SENSE]
const secondFlexibility = collect(0.05:0.05:0.2)
const parsimonious = true

# Boolean parameters
const booleanSpeciesPath = "../../../ModelFiles/Boolean/species.txt"
const booleanRulesPath = "../../../ModelFiles/Boolean/rules.txt"
const TFpath = "../../../ModelFiles/Boolean/TFtargets.txt"

# ODE parameters
const formationRate0 = 0.0001
const repairRate0 = 0.0005
const maxGrowthRate = 0.35

# set initial conditions for ODE
# assume 46% of the drymass are functional proteins (Famili et al., PNAS 2003), 
# there is no damage and the initial drymass is 1 since we are only interested
# in its relativ change, not the absolut 
const P = 0.46
const D = 0.0
const mass = 1.0
const timestep = 0.1
const maxTimesteps = 1000

# parameters for division
const sizeProportion = 0.64
const retention = 0.3

# signalling parameters
const nDelay = 5
const glucoseThreshold = 3.2914
const damageThreshold = 0.001
const trxThreshold = 0.000000002

# regulation parameters
const regulationFactor = collect(0.0:0.02:0.04)

