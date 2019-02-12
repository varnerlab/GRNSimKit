## Solution routines

__GRNSteadyStateSolve__: To estimate the steady-state solution for a GRN system we use the ``GRNSteadyStateSolve`` method.
The ``GRNSteadyStateSolve`` takes a data dictionary of type ``Dict{Symbol,Any}`` and returns the steady-state concentration
vector in your concentration units:

```.jl
    function GRNSteadyStateSolve(data_dictionary::Dict{Symbol,Any})
```
