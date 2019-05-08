
# --- PUBLIC METHODS ----------------------------------------------------------- #
function build_discrete_dynamic_data_dictionary_from_model_file(time_step::Float64, path_to_model_file::String)

    # TODO: Do we have a legit path to the model file?
    check_file_existence(path_to_model_file)

    # switch on the model file type -
    extension = extract_extension(path_to_model_file)
    if (extension == nothing)
        throw(ArgumentError("model file has no extension. GRNSimKit supports models encoded in .toml or .json files"))
    end

    # check the extension -
    extension_symbol = Symbol(extension[2:end])
    if extension_symbol == :toml

        # Load the model dictionary -
        model_dictionary = TOML.parsefile(path_to_model_file)

        # buld the data dictionary and return -
        return build_toml_discrete_dynamic_data_dictionary(time_step, model_dictionary)

    elseif extension_symbol == :json

        # Load the model dictionary -
        model_dictionary = JSON.parsefile(path_to_model_file)

        # buld the data dictionary and return -
        return build_json_discrete_dynamic_data_dictionary(time_step, model_dictionary)

    else
        throw(ArgumentError("unsupported file extension. GRNSimKit supports models encoded in .toml or .json files"))
    end
end

function build_dynamic_data_dictionary_from_model_file(path_to_model_file::String)

    # TODO: Do we have a legit path to the model file?
    check_file_existence(path_to_model_file)

    # switch on the model file type -
    extension = extract_extension(path_to_model_file)
    if (extension == nothing)
        throw(ArgumentError("model file has no extension. GRNSimKit supports models encoded in .toml or .json files"))
    end

    # check the extension -
    extension_symbol = Symbol(extension[2:end])
    if extension_symbol == :toml

        # Load the model dictionary -
        model_dictionary = TOML.parsefile(path_to_model_file)

        # buld the data dictionary and return -
        return build_toml_dynamic_data_dictionary(model_dictionary)

    elseif extension_symbol == :json

        # Load the model dictionary -
        model_dictionary = JSON.parsefile(path_to_model_file)

        # buld the data dictionary and return -
        return build_json_dynamic_data_dictionary(model_dictionary)

    else
        throw(ArgumentError("unsupported file extension. GRNSimKit supports models encoded in .toml or .json files"))
    end
end
# --- PUBLIC METHODS ----------------------------------------------------------- #

# --- INTERNAL METHODS --------------------------------------------------------- #
function build_biophysical_dictionary(model_dictionary::Dict{String,Any})

    # initialize -
    biophysical_dictionary = Dict{Symbol,Any}()

    # convert -
    list_of_constant_dictionaries = model_dictionary["biophysical_constants"]
    for (key, local_dictionary) in list_of_constant_dictionaries

        # create a symbol from the key -
        key_symbol = Symbol(key)

        # get the value -
        value = parse(Float64, local_dictionary["value"])

        # cache -
        biophysical_dictionary[key_symbol] = value
    end

    # return -
    return biophysical_dictionary
end

function build_toml_discrete_dynamic_data_dictionary(time_step::Float64, model_dictionary::Dict{String,Any})
end

function build_toml_dynamic_data_dictionary(time_step::Float64, model_dictionary::Dict{String,Any})
end

function build_json_discrete_dynamic_data_dictionary(time_step::Float64, model_dictionary::Dict{String,Any})

    # build the build_biophysical_dictionary -
    biophysical_dictionary = build_biophysical_dictionary(model_dictionary)

    # build the dilution matrix -
    AM = build_discrete_dilution_matrix(model_dictionary)

    # Build the stoichiometric_matrix -
    SM = build_stoichiometic_matrix(model_dictionary)

    # Compute AHAT and SHAT -
    (nr,nc) = size(AM)
    IM = Matrix(1.0I,nr,nc)
    AHAT = exp(AM*time_step)
    SHAT = inv(AM)*(AHAT - IM)*SM

    # for steady-state calculations -
    CHAT = inv((IM - AHAT))

    # need to correct the AHAT -
    AHAT = correct_discrete_dilution_matrix(AHAT, model_dictionary)

    # setup the system size -
    number_of_states = compute_total_number_of_states(model_dictionary)

    # setup initial conditions -
    initial_condition_array = zeros(number_of_states)

    # compute the transcription rates for all the genes -
    transcription_kinetics_array = compute_transcription_kinetic_array(biophysical_dictionary, model_dictionary)

    # precompute some translation parameters -
    translation_parameters_array = compute_translation_parameter_array(biophysical_dictionary, model_dictionary)

    # precompute a mapping between the species symbols, and index in the x-vector
    species_symbol_index_map = compute_symbol_index_map(model_dictionary)

    # we need bounds for steady state calculations -
    lower_bound_array = zeros(number_of_states)
    upper_bound_array = Inf*ones(number_of_states)

    # --- populate the DD -------------------------------------- #
    data_dictionary = Dict{Symbol,Any}()
    data_dictionary[:problem_type_flag] = :discrete_dynamic
    data_dictionary[:initial_condition_array] = initial_condition_array
    data_dictionary[:number_of_states] = number_of_states
    data_dictionary[:dilution_matrix] = AHAT
    data_dictionary[:stoichiometric_matrix] = SHAT
    data_dictionary[:steady_state_mass_matrix] = CHAT
    data_dictionary[:transcription_kinetics_array] = transcription_kinetics_array
    data_dictionary[:translation_parameters_array] = translation_parameters_array
    data_dictionary[:species_symbol_index_map] = species_symbol_index_map
    data_dictionary[:model_dictionary] = model_dictionary

    # capture the time step -
    data_dictionary[:time_step_size] = time_step

    # return the dd w/default values -
    return data_dictionary
    # ---------------------------------------------------------- #
end

function build_json_dynamic_data_dictionary(model_dictionary::Dict{String,Any})

    # initailzie -
    data_dictionary = Dict{Symbol,Any}()

    # build the build_biophysical_dictionary -
    biophysical_dictionary = build_biophysical_dictionary(model_dictionary)

    # build the dilution matrix -
    AM = build_dilution_matrix(model_dictionary)

    # Build the stoichiometric_matrix -
    SM = build_stoichiometic_matrix(model_dictionary)

    # setup the system size -
    number_of_states = compute_total_number_of_states(model_dictionary)

    # setup initial conditions -
    initial_condition_array = zeros(number_of_states)

    # compute the transcription rates for all the genes -
    transcription_kinetics_array = compute_transcription_kinetic_array(biophysical_dictionary, model_dictionary)

    # precompute some translation parameters -
    translation_parameters_array = compute_translation_parameter_array(biophysical_dictionary, model_dictionary)

    # precompute a mapping between the species symbols, and index in the x-vector
    species_symbol_index_map = compute_symbol_index_map(model_dictionary)

    # --- populate the DD -------------------------------------- #
    data_dictionary[:problem_type_flag] = :general_dynamic
    data_dictionary[:initial_condition_array] = initial_condition_array
    data_dictionary[:number_of_states] = number_of_states
    data_dictionary[:dilution_matrix] = AM
    data_dictionary[:stoichiometric_matrix] = SM
    data_dictionary[:transcription_kinetics_array] = transcription_kinetics_array
    data_dictionary[:translation_parameters_array] = translation_parameters_array
    data_dictionary[:species_symbol_index_map] = species_symbol_index_map
    data_dictionary[:model_dictionary] = model_dictionary

    # return the dd w/default values -
    return data_dictionary
end

function build_differential_algebraic_data_dictionary(path_to_model_file::String)

    # initailzie -
    data_dictionary = Dict{Symbol,Any}()

    # Load the model dictionary -
    model_dictionary = JSON.parsefile(path_to_model_file)

    # build the build_biophysical_dictionary -
    biophysical_dictionary = build_biophysical_dictionary(model_dictionary)

    # compute the transcription rates for all the genes -
    transcription_kinetics_array = compute_transcription_kinetic_array(biophysical_dictionary, model_dictionary)

    # precompute some translation parameters -
    translation_parameters_array = compute_translation_parameter_array(biophysical_dictionary, model_dictionary)

    # precompute a mapping between the species symbols, and index in the x-vector
    species_symbol_index_map = compute_symbol_index_map(model_dictionary)

    # build the dilution matrix -
    AM = build_dilution_matrix(model_dictionary)

    # Build the stoichiometric_matrix -
    SM = build_stoichiometic_matrix(model_dictionary)

    # --- populate the DD ------------------------------------------------ #
    data_dictionary[:problem_type_flag] = :differential_algebraic_dynamic
    data_dictionary[:transcription_kinetics_array] = transcription_kinetics_array
    data_dictionary[:translation_parameters_array] = translation_parameters_array
    data_dictionary[:dilution_matrix] = AM
    data_dictionary[:stoichiometric_matrix] = SM

    data_dictionary[:species_symbol_index_map] = species_symbol_index_map
    data_dictionary[:model_dictionary] = model_dictionary

    # return the dd w/default values -
    return data_dictionary
    # ------------------------------------------------------------------- #
end
# --- INTERNAL METHODS --------------------------------------------------------- #

# --- DEPRECATED ------------------------------------------------------------------------ #
function build_default_data_dictionary(path_to_model_file::String)

    # TODO: Do we have a legit path to the model file?
    # ....

    # initailzie -
    data_dictionary = Dict{Symbol,Any}()

    # Load the model dictionary -
    model_dictionary = JSON.parsefile(path_to_model_file)

    # build the build_biophysical_dictionary -
    biophysical_dictionary = build_biophysical_dictionary(model_dictionary)

    # build the dilution matrix -
    AM = build_dilution_matrix(model_dictionary)

    # Build the stoichiometric_matrix -
    SM = build_stoichiometic_matrix(model_dictionary)

    # setup the system size -
    number_of_states = compute_total_number_of_states(model_dictionary)

    # setup initial conditions -
    initial_condition_array = zeros(number_of_states)

    # compute the transcription rates for all the genes -
    transcription_kinetics_array = compute_transcription_kinetic_array(biophysical_dictionary, model_dictionary)

    # precompute some translation parameters -
    translation_parameters_array = compute_translation_parameter_array(biophysical_dictionary, model_dictionary)

    # precompute a mapping between the species symbols, and index in the x-vector
    species_symbol_index_map = compute_symbol_index_map(model_dictionary)

    # --- populate the DD -------------------------------------- #
    data_dictionary[:problem_type_flag] = :general_dynamic
    data_dictionary[:initial_condition_array] = initial_condition_array
    data_dictionary[:number_of_states] = number_of_states
    data_dictionary[:dilution_matrix] = AM
    data_dictionary[:stoichiometric_matrix] = SM
    data_dictionary[:transcription_kinetics_array] = transcription_kinetics_array
    data_dictionary[:translation_parameters_array] = translation_parameters_array
    data_dictionary[:species_symbol_index_map] = species_symbol_index_map
    data_dictionary[:model_dictionary] = model_dictionary

    # return the dd w/default values -
    return data_dictionary
    # ---------------------------------------------------------- #
end
# --- DEPRECATED ------------------------------------------------------------------------ #
