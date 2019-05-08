module GRNSimKit

# include -
include("Include.jl")

# exports -
export GRNSteadyStateSolve
export GRNDynamicSolve
export GRNDiscreteDynamicSolve
export GRNDifferentialAlgebraicSolve

export build_default_data_dictionary
export build_discrete_dynamic_data_dictionary
export build_discrete_dynamic_data_dictionary_from_model_file
export build_dynamic_data_dictionary_from_model_file

export write_discrete_perturbed_simulation_files
export write_scaled_sensitivity_array
export time_average_senstivity_array

end # module
