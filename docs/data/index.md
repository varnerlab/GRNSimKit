## Data dictionary
Model parameters, and model structure are wrapped into a data dictionary that is used by the simulation routines
to solve the model equations. You create a data dictionary by calling one of the ``build_{*}_data_dictionary`` methods
(where ``{*}`` is replaced by a simulation type) with the path to the [model specification file](../model/index.md) as an argument.  

### Functions
[build_default_data_dictionary](https://github.com/varnerlab/GRNSimKit/blob/master/src/Data.jl): Takes the path to a
model specification file (String) and input, and returns a Julia dictionary (Dict{Symbol,Any}) with parameter and model structure information.

```jl

  function build_default_data_dictionary(path_to_model_file::String) -> Dict{Symbol,Any}

``` 
