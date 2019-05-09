function extract_extension(url::String)::Union{String,Nothing}
    try
        return match(r"\.[A-Za-z0-9]+$", url).match
    catch
        return nothing
    end
end

function check_file_existence(url::String)

    if (isfile(url) == false)
        throw(ArgumentError("opening file $(url): No such file or directory"))
    end
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
    number_of_nodes = length(list_of_gene_nodes)
    for node_index = 1:number_of_nodes

        # get the dictionary for this node -
        node_key = "N$(node_index)"
        node_dictionary = list_of_gene_nodes[node_key]

        # grab the parameters for this node -
        mRNA_half_life_in_hr = parse(Float64, node_dictionary["parameters"]["mRNA_half_life_in_hr"])

        # compute the degradation rate constant -
        kDX = -(1/mRNA_half_life_in_hr)*log(0.5)

        # dilution term -
        term = (mu+kDX)

        # push -
        push!(dilution_term_array, term)
    end

    for node_index = 1:number_of_nodes

        # get the dictionary for this node -
        node_key = "N$(node_index)"
        node_dictionary = list_of_gene_nodes[node_key]

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

function build_toml_discrete_dilution_matrix(model_dictionary::Dict{String,Any})

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
    tmp_value = model_dictionary["biophysical_constants"]["cell_doubling_time"]
    cell_doubling_time = parse(Float64, tmp_value)

    # compute the mumax -
    mu = log(2)/cell_doubling_time

    # compute the rate array -
    dilution_term_array = Float64[]

    # ok, next let's compute the degrdation for each mRNA and protein node -
    list_of_gene_nodes = model_dictionary["list_of_nodes"]
    number_of_nodes = length(list_of_gene_nodes)
    for node_index = 1:number_of_nodes

        # get the dictionary for this node -
        node_key = "N$(node_index)"
        node_dictionary = list_of_gene_nodes[node_key]

        # grab the parameters for this node -
        mRNA_half_life_in_hr = parse(Float64, node_dictionary["parameters"]["mRNA_half_life_in_hr"])

        # compute the degradation rate constant -
        kDX = -(1/mRNA_half_life_in_hr)*log(0.5)

        # dilution term -
        term = (mu+kDX)

        # push -
        push!(dilution_term_array, term)
    end

    for node_index = 1:number_of_nodes

        # get the dictionary for this node -
        node_key = "N$(node_index)"
        node_dictionary = list_of_gene_nodes[node_key]

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

function build_toml_dilution_matrix(model_dictionary::Dict{String,Any})

    # get node and species lists -
    list_of_species_dictionaries = model_dictionary["list_of_species_symbols"]

    # how many species do we have?
    total_number_of_species = length(list_of_species_dictionaries)

    # build a blank st array -
    dilution_matrix = -1*Matrix{Float64}(I, total_number_of_species, total_number_of_species)

    # correct -
    for species_dictionary in list_of_species_dictionaries

        # get the index, and the type -
        species_index = species_dictionary["index"]
        species_type_symbol = Symbol(species_dictionary["type"])

        # if symbol is of type constant, then the dilution element is 0
        if species_type_symbol == :constant
            dilution_matrix[species_index,species_index] = 0.0
        end
    end

    # get the cell doubling time -
    tmp_value = model_dictionary["biophysical_constants"]["cell_doubling_time"]
    cell_doubling_time = parse(Float64, tmp_value)

    # compute the mumax -
    mu = log(2)/cell_doubling_time

    # compute the rate array -
    dilution_term_array = Float64[]

    # ok, next let's compute the degrdation for each mRNA and protein node -
    list_of_gene_nodes = model_dictionary["list_of_nodes"]
    number_of_nodes = length(list_of_gene_nodes)
    for node_index = 1:number_of_nodes

        # get the dictionary for this node -
        node_key = "N$(node_index)"
        node_dictionary = list_of_gene_nodes[node_key]

        # grab the parameters for this node -
        mRNA_half_life_in_hr = parse(Float64, node_dictionary["parameters"]["mRNA_half_life_in_hr"])

        # compute the degradation rate constant -
        kDX = -(1/mRNA_half_life_in_hr)*log(0.5)

        # dilution term -
        term = (mu+kDX)

        # push -
        push!(dilution_term_array, term)
    end

    for node_index = 1:number_of_nodes

        # get the dictionary for this node -
        node_key = "N$(node_index)"
        node_dictionary = list_of_gene_nodes[node_key]

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
    number_of_nodes = length(list_of_gene_nodes)
    for node_index = 1:number_of_nodes

        # get the dictionary for this node -
        node_key = "N$(node_index)"
        node_dictionary = list_of_gene_nodes[node_key]

        # grab the parameters for this node -
        mRNA_half_life_in_hr = parse(Float64, node_dictionary["parameters"]["mRNA_half_life_in_hr"])

        # compute the degradation rate constant -
        kDX = -(1/mRNA_half_life_in_hr)*log(0.5)

        # dilution term -
        term = (mu+kDX)

        # push -
        push!(dilution_term_array, term)
    end

    for node_index = 1:number_of_nodes

        # get the dictionary for this node -
        node_key = "N$(node_index)"
        node_dictionary = list_of_gene_nodes[node_key]

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

function compute_total_number_of_states(model_dictionary::Dict{String,Any})

    # get node and species lists -
    list_of_species_dictionaries = model_dictionary["list_of_species_symbols"]

    # how many species do we have?
    total_number_of_species = length(list_of_species_dictionaries)

    # return -
    return total_number_of_species
end
