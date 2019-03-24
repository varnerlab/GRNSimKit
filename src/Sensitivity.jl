function get_parameter_value(parameter_path_array::Array{String,1}, model_dictionary::Dict{String, Any})

    # get a key symbol -
    next_key_symbol = pop!(parameter_path_array)

    # recursive descend -
    if (typeof(get(model_dictionary,next_key_symbol,Nothing)) == Dict{String,Any})

        # get dictionary, and go down again ...
        next_dictionary = model_dictionary[next_key_symbol]
        return get_parameter_value(parameter_path_array, next_dictionary)
    elseif (typeof(get(model_dictionary,next_key_symbol,Nothing)) == String)

        # return the value -
        return model_dictionary[next_key_symbol]
    else
        return Nothing
    end
end

function set_parameter_value(parameter_path_array::Array{String,1}, model_dictionary::Dict{String, Any}, value::String)

    # get a key symbol -
    next_key_symbol = pop!(parameter_path_array)

    # recursive descend -
    if (typeof(get(model_dictionary,next_key_symbol,Nothing)) == Dict{String,Any})

        # get dictionary, and go down again ...
        next_dictionary = model_dictionary[next_key_symbol]
        return set_parameter_value(parameter_path_array, next_dictionary, value)

    elseif (typeof(get(model_dictionary,next_key_symbol,Nothing)) == String)

        # set the new value -
        model_dictionary[next_key_symbol] = value
        return true
    else
        return false
    end
end

function perturb_model_dictionary(parameter_path::String, model_dictionary::Dict{String,Any}, perturbation_direction::Symbol)

    # create a copy the model dictionary -
    copy_model_dictionary = deepcopy(model_dictionary)

    # get the value -
    parameter_path_array = reverse(splitpath(parameter_path))
    nominal_value = get_parameter_value(parameter_path_array, model_dictionary)
    if nominal_value == Nothing
        error_msg = "No value at $(parameter_path). Please check the path"
        error(error_msg)
    end

    # ok, if we get here we have a value -
    if perturbation_direction == :UP

        # crete a new value for this parameter -
        nominal_value_float = parse(Float64,nominal_value)
        new_value = nominal_value_float*(1+0.05)

        # set new value of the dictionary copy -
        parameter_path_array = reverse(splitpath(parameter_path))
        did_set = set_parameter_value(parameter_path_array, copy_model_dictionary, string(new_value))
        if (did_set == false)
            error_msg = "Failed to set new value at $(parameter_path)"
            error(error_msg)
        end

        # return -
        return copy_model_dictionary

    elseif perturbation_direction == :DOWN

        # crete a new value for this parameter -
        nominal_value_float = parse(Float64,nominal_value)
        new_value = nominal_value_float*(1-0.05)

        # set new value of the dictionary copy -
        parameter_path_array = reverse(splitpath(parameter_path))
        did_set = set_parameter_value(parameter_path_array, copy_model_dictionary, string(new_value))
        if (did_set == false)
            error_msg = "Failed to set new value at $(parameter_path)"
            error(error_msg)
        end

        # return -
        return copy_model_dictionary
    else
        error_msg = "Perturbation parameter is either :UP or :DOWN. You asked for :$(perturbation_direction)"
        error(error_msg)
    end
end

function write_discrete_perturbed_simulation_files(time_span::Tuple{Float64, Float64, Float64}, initial_condition_array::Array{Float64,1}, parameter_path_array::Array{String}, model_dictionary::Dict{String,Any}, output_dir::String)

    # get the time_step_size
    time_step_size = time_span[end]

    # build_discrete_dynamic_data_dictionary -
    nominal_dictionary = build_discrete_dynamic_data_dictionary(time_step_size, model_dictionary)

    # set the nominal steady-state as the initial condition -
    nominal_dictionary[:initial_condition_array] = initial_condition_array

    # solve the model in the nominal case -
    (TN, XN) = GRNDiscreteDynamicSolve(time_span, nominal_dictionary)

    # dump nominal simulation to disk -
    nominal_simulation = [TN XN]
    output_file_path = joinpath(output_dir,"simulation_nominal.dat")
    writedlm(output_file_path, nominal_simulation)

    # ok, so now lets do the perturbations -
    counter = 1
    for parameter_path in parameter_path_array

        # create perturbed dictionaries -
        perturbed_up_model_dictionary = perturb_model_dictionary(parameter_path, model_dictionary, :UP)
        perturbed_down_model_dictionary = perturb_model_dictionary(parameter_path, model_dictionary, :DOWN)

        # ok, evaulate the up case -
        up_dictionary = build_discrete_dynamic_data_dictionary(time_step_size, perturbed_up_model_dictionary)
        up_dictionary[:initial_condition_array] = initial_condition_array
        (TU, XU) = GRNDiscreteDynamicSolve(time_span, up_dictionary)
        file_name = "simulation_up_P$(counter).dat"
        output_file_path = joinpath(output_dir,file_name)
        writedlm(output_file_path, [TU XU])

        # ok, evaulate the down case -
        down_dictionary = build_discrete_dynamic_data_dictionary(time_step_size, perturbed_down_model_dictionary)
        down_dictionary[:initial_condition_array] = initial_condition_array
        (TD, XD) = GRNDiscreteDynamicSolve(time_span, up_dictionary)
        file_name = "simulation_down_P$(counter).dat"
        output_file_path = joinpath(output_dir,file_name)
        writedlm(output_file_path, [TD XD])

        # update the counter -
        counter = counter + 1
    end
end
