## Solution routines

__GRNSteadyStateSolve__: To estimate the steady-state solution for a GRN system we use the ``GRNSteadyStateSolve`` method:

```jl
    GRNSteadyStateSolve: Estimate the steady-state concentration vector for a gene regulatory network
    input arguments: data_dictionary::Dict{Symbol,Any}
    output arguments: X::Array{Float64,1}
```
