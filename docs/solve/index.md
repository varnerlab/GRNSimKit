## GRNSimKit solution routines
There are several ways to solve GRN models, depending upon what you want to do.

* [GRNSteadyStateSolve](https://github.com/varnerlab/GRNSimKit/blob/master/src/Solve.jl): Estimates the steady-state solution for a GRN system. The ``GRNSteadyStateSolve`` method takes a data dictionary [Dict{Symbol,Any}](https://docs.julialang.org/en/v1/base/collections/#Dictionaries-1) holding system parameters, and returns the steady-state concentration
vector (in your concentration units) [Array{Float64,1}](https://docs.julialang.org/en/v1/base/arrays/):

```jl
    function GRNSteadyStateSolve(data_dictionary::Dict{Symbol,Any}) -> Array{Float64,1}
```

* [GRNDynamicSolve](https://github.com/varnerlab/GRNSimKit/blob/master/src/Solve.jl):
Solves the system of ordinary differential equations using the routines in the [DifferentialEquations](https://github.com/JuliaDiffEq/DifferentialEquations.jl) package.

```jl
  function GRNDynamicSolve(time_span::Tuple{Float64,Float64}, data_dictionary::Dict{Symbol,Any}) -> (T::Array{Float64,1}, X::Array{Float64,2})
```
