function compute_total_number_of_states(model_dictionary::Dict{String,Any})

    # get node and species lists -
    list_of_species_dictionaries = model_dictionary["list_of_species_symbols"]

    # how many species do we have?
    total_number_of_species = length(list_of_species_dictionaries)

    # return -
    return total_number_of_species
end

function compute_transcription_kinetic_array(biophysical_dictionary::Dict{Symbol,Any}, model_dictionary::Dict{String,Any})

    # get stuff from the biophysical dictionary -
    RNAPII_copy_number = biophysical_dictionary[:copies_of_rnapII_per_cell]
    characteristic_length = biophysical_dictionary[:characteristic_transcript_length]
    transcription_elongation_rate = biophysical_dictionary[:transcription_elongation_rate]
    characteristic_initiation_time = biophysical_dictionary[:characteristic_initiation_time_transcription]
    KX = biophysical_dictionary[:transcription_saturation_constant]
    mass_of_single_cell = biophysical_dictionary[:mass_of_single_cell]
    fraction_of_water_per_cell = biophysical_dictionary[:fraction_of_water_per_cell]

    # fraction of dry weight -
    fraction_dry_cell = 1 - fraction_of_water_per_cell      # dimensionless

    # avagodros number -
    av_number = 6.02e23                                     # number/mol

    # what is the RNAP concentration -
    RNAPII_concentration = RNAPII_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell) # nmol/gDW

    # initialize -
    transcription_kinetics_array = Float64[]

    # compute the kinetic limit -
    list_of_gene_nodes = model_dictionary["list_of_nodes"]
    for node_dictionary in list_of_gene_nodes

        # compute elongation constant -
        gene_length = parse(Float64,node_dictionary["parameters"]["gene_length_in_nt"])
        gene_copy_number = parse(Float64, node_dictionary["parameters"]["gene_copy_number_per_cell"])

        # compute kE -
        kE = transcription_elongation_rate*(1/gene_length)

        # compute kI -
        kI = (1/characteristic_initiation_time)

        # compute the tau factor -
        tau_factor = (kE/kI)

        # compute the gene concentration -
        G = gene_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell)   # nmol/gDW

        # Compute the rate -
        value = kE*RNAPII_concentration*(G/(KX*tau_factor+(1+tau_factor)*G))*(3600)

        # push -
        push!(transcription_kinetics_array, value)
    end

    # return -
    return transcription_kinetics_array
end

function compute_translation_parameter_array(biophysical_dictionary::Dict{Symbol,Any}, model_dictionary::Dict{String,Any})

    # get some stuff from the biophysical_dictionary -
    mass_of_single_cell = biophysical_dictionary[:mass_of_single_cell]
    fraction_of_water_per_cell = biophysical_dictionary[:fraction_of_water_per_cell]
    ribosome_copy_number = biophysical_dictionary[:copies_of_ribosome_per_cell]
    translation_elongation_rate = biophysical_dictionary[:translation_elongation_rate]
    characteristic_initiation_time_translation = biophysical_dictionary[:characteristic_initiation_time_translation]
    KL = biophysical_dictionary[:translation_saturation_constant]

    # fraction of dry weight -
    fraction_dry_cell = 1 - fraction_of_water_per_cell      # dimensionless

    # avagodros number -
    av_number = 6.02e23                                     # number/mol

    # what is the ribosome concentration -
    ribosome_concentration = ribosome_copy_number*(1/mass_of_single_cell)*(1/av_number)*(1e9)*(1/fraction_dry_cell) # nmol/gDW

    # initialize -
    translation_parameters_array = Array{TranslationParameters,1}()

    # compute the translation parameters -
    list_of_gene_nodes = model_dictionary["list_of_nodes"]
    for node_dictionary in list_of_gene_nodes

        # compute elongation constant -
        protein_length = parse(Float64,node_dictionary["parameters"]["protein_length_in_aa"])

        # compute kE -
        kE = translation_elongation_rate*(1/protein_length)

        # compute kI -
        kI = (1/characteristic_initiation_time_translation)

        # compute the tau factor -
        tau_factor = (kE/kI)

        # compute the vmax -
        translation_vmax = kE*ribosome_concentration

        # build -
        translationParameters = TranslationParameters()
        translationParameters.vmax = translation_vmax*(3600)    # convert to hr
        translationParameters.tau_factor = tau_factor
        translationParameters.KL = KL
        push!(translation_parameters_array, translationParameters)
    end

    return translation_parameters_array
end

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

function correct_discrete_dilution_matrix(AHAT,model_dictionary::Dict{String,Any})

    # get node and species lists -
    list_of_species_dictionaries = model_dictionary["list_of_species_symbols"]

    # how many species do we have?
    total_number_of_species = length(list_of_species_dictionaries)

    # build a blank st array -
    dilution_matrix = -1*Matrix{Float64}(I, total_number_of_species, total_number_of_species)

    # correct -
    for species_dictionary in list_of_species_dictionaries

        # get the index, and the type -
        species_index = parse(Int64, species_dictionary["index"])
        species_type_symbol = Symbol(species_dictionary["type"])

        # if symbol is of type constant, then the dilution element is 0
        if species_type_symbol == :constant
            AHAT[species_index,species_index] = 1.0
        end
    end

    return AHAT
end

function build_discrete_dilution_matrix(model_dictionary::Dict{String,Any})

    # get node and species lists -
    list_of_species_dictionaries = model_dictionary["list_of_species_symbols"]

    # how many species do we have?
    total_number_of_species = length(list_of_species_dictionaries)

    # build a blank st array -
    dilution_matrix = -1*Matrix{Float64}(I, total_number_of_species, total_number_of_species)

    # correct -
    for species_dictionary in list_of_species_dictionaries

        # get the index, and the type -
        species_index = parse(Int64, species_dictionary["index"])
        species_type_symbol = Symbol(species_dictionary["type"])

        # if symbol is of type constant, then the dilution element is 0
        if species_type_symbol == :constant
            dilution_matrix[species_index,species_index] = 1.0
        end
    end

    # get the cell doubling time -
    tmp_value = model_dictionary["biophysical_constants"]["cell_doubling_time"]["value"]
    cell_doubling_time = parse(Float64, tmp_value)

    # compute the mumax -
    mu = log(2)/cell_doubling_time

    # compute the rate array -
    dilution_term_array = Float64[]

    # ok, next let's compute the degrdation for each mRNA and protein node -
    list_of_gene_nodes = model_dictionary["list_of_nodes"]
    for node_dictionary in list_of_gene_nodes

        # grab the parameters for this node -
        mRNA_half_life_in_hr = parse(Float64, node_dictionary["parameters"]["mRNA_half_life_in_hr"])

        # compute the degradation rate constant -
        kDX = -(1/mRNA_half_life_in_hr)*log(0.5)

        # dilution term -
        term = (mu+kDX)

        # push -
        push!(dilution_term_array, term)
    end

    for node_dictionary in list_of_gene_nodes

        # grab the parameters for this node -
        protein_half_life_in_hr = parse(Float64, node_dictionary["parameters"]["protein_half_life_in_hr"])

        # compute the degradation rate constant -
        kDL = -(1/protein_half_life_in_hr)*log(0.5)

        # dilution term -
        term = (mu+kDL)

        # push -
        push!(dilution_term_array, term)
    end

    # add zeros for constant terms -
    for species_dictionary in list_of_species_dictionaries

        # get the type -
        species_type_symbol = Symbol(species_dictionary["type"])

        # if symbol is of type constant, then the dilution element is 0
        if species_type_symbol == :constant
            push!(dilution_term_array, 1.0)
        end
    end

    # return the dilution array -
    return dilution_matrix.*dilution_term_array
end

function build_dilution_matrix(model_dictionary::Dict{String,Any})

    # get node and species lists -
    list_of_species_dictionaries = model_dictionary["list_of_species_symbols"]

    # how many species do we have?
    total_number_of_species = length(list_of_species_dictionaries)

    # build a blank st array -
    dilution_matrix = -1*Matrix{Float64}(I, total_number_of_species, total_number_of_species)

    # correct -
    for species_dictionary in list_of_species_dictionaries

        # get the index, and the type -
        species_index = parse(Int64, species_dictionary["index"])
        species_type_symbol = Symbol(species_dictionary["type"])

        # if symbol is of type constant, then the dilution element is 0
        if species_type_symbol == :constant
            dilution_matrix[species_index,species_index] = 0.0
        end
    end

    # get the cell doubling time -
    tmp_value = model_dictionary["biophysical_constants"]["cell_doubling_time"]["value"]
    cell_doubling_time = parse(Float64, tmp_value)

    # compute the mumax -
    mu = log(2)/cell_doubling_time

    # compute the rate array -
    dilution_term_array = Float64[]

    # ok, next let's compute the degrdation for each mRNA and protein node -
    list_of_gene_nodes = model_dictionary["list_of_nodes"]
    for node_dictionary in list_of_gene_nodes

        # grab the parameters for this node -
        mRNA_half_life_in_hr = parse(Float64, node_dictionary["parameters"]["mRNA_half_life_in_hr"])

        # compute the degradation rate constant -
        kDX = -(1/mRNA_half_life_in_hr)*log(0.5)

        # dilution term -
        term = (mu+kDX)

        # push -
        push!(dilution_term_array, term)
    end

    for node_dictionary in list_of_gene_nodes

        # grab the parameters for this node -
        protein_half_life_in_hr = parse(Float64, node_dictionary["parameters"]["protein_half_life_in_hr"])

        # compute the degradation rate constant -
        kDL = -(1/protein_half_life_in_hr)*log(0.5)

        # dilution term -
        term = (mu+kDL)

        # push -
        push!(dilution_term_array, term)
    end

    # add zeros for constant terms -
    for species_dictionary in list_of_species_dictionaries

        # get the type -
        species_type_symbol = Symbol(species_dictionary["type"])

        # if symbol is of type constant, then the dilution element is 0
        if species_type_symbol == :constant
            push!(dilution_term_array, 0.0)
        end
    end

    # return the dilution array -
    return dilution_matrix.*dilution_term_array
end

function build_stoichiometic_matrix(model_dictionary::Dict{String,Any})

    # get node and species lists -
    list_of_species_dictionaries = model_dictionary["list_of_species_symbols"]
    list_of_gene_nodes = model_dictionary["list_of_nodes"]

    # how many species do we have?
    number_of_species = length(list_of_species_dictionaries)
    number_of_reactions = 2*length(list_of_gene_nodes)

    # build a blank st array -
    stoichiometric_matrix = zeros(number_of_species, number_of_reactions)

    # fill in the missing 1's -
    for species_index = 1:number_of_species
        for reaction_index = 1:number_of_reactions

            if (species_index == reaction_index)
                stoichiometric_matrix[species_index,reaction_index] = 1.0
            end
        end
    end

    # return -
    return stoichiometric_matrix
end

function compute_symbol_index_map(model_dictionary)

    # initlaize -
    species_index_map = Dict{String,Int}()

    # get the list of species -
    list_of_species_symbols = model_dictionary["list_of_species_symbols"]
    for species_object_dict in list_of_species_symbols

        # get the symbol -
        species_symbol = species_object_dict["symbol"]
        species_index = parse(Int64, species_object_dict["index"])

        # cache -
        species_index_map[species_symbol] = species_index
    end

    # return -
    return species_index_map
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

function build_discrete_dynamic_data_dictionary(time_step::Float64, model_dictionary::Dict{String,Any})

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

function build_discrete_dynamic_data_dictionary_from_model_file(time_step::Float64, path_to_model_file::String)

    # TODO: Do we have a legit path to the model file?
    # ....

    # Load the model dictionary -
    model_dictionary = JSON.parsefile(path_to_model_file)

    # buld the data dictionary and return -
    return build_discrete_dynamic_data_dictionary(time_step, model_dictionary)
end

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
