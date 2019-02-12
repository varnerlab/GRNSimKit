## Solution routines

__Steady state solution__: To estimate the steady-state solution for a GRN system we use the ``GRNSteadyStateSolve`` method:

     GRNSteadyStateSolve:
     input arguments: data_dictionary::Dict{Symbol,Any}
     output arguments: X::Array{Float64,1}
