
function calculate_control_array(t,x,data_dictionary)

    # initialize -
    control_array = Float64[]

    # grab the model_dictionary -
    model_dictionary = data_dictionary[:model_dictionary]

    # grab the species mapping -
    species_symbol_index_map = data_dictionary[:species_symbol_index_map]

    # grab the list of nodes -
    list_of_nodes = model_dictionary["list_of_nodes"]
    for node in list_of_nodes

        # initialize -
        numerator_term_array = Float64[]
        denominator_term_array = Float64[]

        # --------------------------------------------------------------------------------- #
        # background -
        W_background = parse(Float64, node["parameters"]["background_expression"])
        push!(numerator_term_array, W_background)

        # list of activators -
        list_of_activators = node["transcription_control_structure"]["list_of_activators"]
        for activator_node in list_of_activators

            # grab the symbol and stuff -
            species_symbol = activator_node["symbol"]
            species_index = species_symbol_index_map[species_symbol]
            state_value = abs(x[species_index])

            # binding stuff -
            W = parse(Float64, activator_node["weight"])
            n = parse(Float64, activator_node["binding_function"]["order"])
            K = parse(Float64, activator_node["binding_function"]["binding_constant"])

            # compute the binding function -
            f_binding = ((state_value/K)^n)/(1+(state_value/K)^n)
            push!(numerator_term_array, W*f_binding)
        end

        # list of repressors -
        list_of_repressors = node["transcription_control_structure"]["list_of_repressors"]
        push!(denominator_term_array, 1.0)
        for repressor_node in list_of_repressors

            # grab the symbol and stuff -
            species_symbol = repressor_node["symbol"]
            species_index = species_symbol_index_map[species_symbol]
            state_value = abs(x[species_index])

            # binding stuff -
            W = parse(Float64, repressor_node["weight"])
            n = parse(Float64, repressor_node["binding_function"]["order"])
            K = parse(Float64, repressor_node["binding_function"]["binding_constant"])

            # compute the binding function -
            f_binding = ((state_value/K)^n)/(1+(state_value/K)^n)
            push!(denominator_term_array, W*f_binding)
        end
        # ------------------------------------------------------------------------- #

        # compute the control value -
        value = sum(numerator_term_array)/(sum(numerator_term_array)+sum(denominator_term_array))

        # cache -
        push!(control_array, value)
    end

    # add 1's for translation -
    for node in list_of_nodes
        push!(control_array, 1.0)
    end

    # return control_array -
    return control_array
end

function calculate_txtl_kinetics_array(t,x,data_dictionary)

    # we *precalculated* the transcription kinetics array -
    transcription_kinetics_array = data_dictionary[:transcription_kinetics_array]

    # we need to calculate the translation rate on the fly
    # since it depends upon mRNA -
    translation_kinetics_array = Float64[]

    # grab the translation parameters from the data dictionary -
    translation_parameters_array = data_dictionary[:translation_parameters_array]
    number_of_translation_rates = length(translation_parameters_array)
    counter = 1
    for index = 1:number_of_translation_rates

        # mRNA is always the first block of species ...
        mRNA_concentration = x[index]

        # grab the parameters struct -
        parameters_struct = translation_parameters_array[index]
        vmax = parameters_struct.vmax
        tau = parameters_struct.tau_factor
        KL = parameters_struct.KL

        # build the rate -
        value = vmax*(mRNA_concentration/(KL*tau+(1+tau)*mRNA_concentration))

        # package -
        push!(translation_kinetics_array, value)
    end

    # finally, build the kinetics array -
    txtl_kinetics_array = Float64[]
    for value in transcription_kinetics_array
        push!(txtl_kinetics_array, value)
    end

    for value in translation_kinetics_array
        push!(txtl_kinetics_array, value)
    end

    # return -
    return txtl_kinetics_array
end

function discrete_steady_state_balance_residual(x,data_dictionary)

    # get the structure arrays from the data_dictionary -
    AHAT = data_dictionary[:dilution_matrix]
    SHAT = data_dictionary[:stoichiometric_matrix]
    CHAT = data_dictionary[:steady_state_mass_matrix]

    # calculate the kinetics array -
    kinetics_array = calculate_txtl_kinetics_array(0,x,data_dictionary)

    # calculate the control array -
    control_array = calculate_control_array(0,x,data_dictionary)

    # modfiy the kinetics -
    rV = kinetics_array.*control_array

    # compute the error -
    residual = (x - CHAT*SHAT*rV)
    error = transpose(residual)*residual

    # return the error
    return error
end

function discrete_dynamic_balances(t,x,data_dictionary)

    # get the structure arrays from the data_dictionary -
    AHAT = data_dictionary[:dilution_matrix]
    SHAT = data_dictionary[:stoichiometric_matrix]

    # calculate the kinetics array -
    kinetics_array = calculate_txtl_kinetics_array(t,x,data_dictionary)

    # calculate the control array -
    control_array = calculate_control_array(t,x,data_dictionary)

    # modfiy the kinetics -
    rV = kinetics_array.*control_array

    # compute the next time step -
    next_step = AHAT*x + SHAT*rV

    # return -
    return next_step
end

function balances(dx,x,data_dictionary,t)

    # get the structure arrays from the data_dictionary -
    AM = data_dictionary[:dilution_matrix]
    SM = data_dictionary[:stoichiometric_matrix]

    # what is my system size?
    number_of_states = data_dictionary[:number_of_states]

    # calculate the kinetics array -
    kinetics_array = calculate_txtl_kinetics_array(t,x,data_dictionary)

    # calculate the control array -
    control_array = calculate_control_array(t,x,data_dictionary)

    # modfiy the kinetics -
    rV = kinetics_array.*control_array

    # compute the dxdt vector -
    dxdt = AM*x+SM*rV

    # package -
    for index = 1:number_of_states
        dx[index] = dxdt[index]
    end
end
