using MAT
using JuMP
using Gurobi
using DelimitedFiles
using Statistics

include("InitialiseJuMPModel.jl")
include("FBA.jl")

include("InitialiseBooleanModel.jl")
include("Boolean.jl")

include("InitialiseODEModel.jl")
include("DamageODE.jl")

include("Integration.jl")
include("Cell.jl")