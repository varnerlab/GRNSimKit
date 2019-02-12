## Solution routines

__GRNSteadyStateSolve__: To estimate the steady-state solution for a GRN system we use the ``GRNSteadyStateSolve`` method.
The ``GRNSteadyStateSolve`` takes a data dictionary of type [Dict{Symbol,Any}](https://docs.julialang.org/en/v1/base/collections/#Dictionaries-1) and returns the steady-state concentration
vector (in your concentration units) of type [Array{Float64,1}](https://docs.julialang.org/en/v1/base/arrays/):

```.jl
    function GRNSteadyStateSolve(data_dictionary::Dict{Symbol,Any})
```
