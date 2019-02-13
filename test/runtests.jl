using GRNSimKit
using Test

# Data tests -
@testset "DefaultDataDictionary" begin

    # define my top level path -
    path_to_package = dirname(pathof(@__MODULE__))

    # path to the model file -
    path_to_model_file = joinpath(path_to_package,"..", "test/models/CT1.json")

    # build the dd -
    dd = build_default_data_dictionary(path_to_model_file)

    # load -
    @test dd[:number_of_states] == 7
end
