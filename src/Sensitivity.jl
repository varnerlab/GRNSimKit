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
        (TD, XD) = GRNDiscreteDynamicSolve(time_span, down_dictionary)
        file_name = "simulation_down_P$(counter).dat"
        output_file_path = joinpath(output_dir,file_name)
        writedlm(output_file_path, [TD XD])

        # update the counter -
        counter = counter + 1
    end
end

function write_scaled_sensitivity_array(time_index_array, state_index_array, parameter_path_array, simulation_data_path, sensitivity_output_dir, model_dictionary::Dict{String,Any})

    # get node and species lists -
    list_of_species_dictionaries = model_dictionary["list_of_species_symbols"]

    # ok, setup some dimensions -
    number_of_parameters = length(parameter_path_array)
    number_of_states = length(state_index_array)
    number_of_time_steps = length(time_index_array)

    # initialize ensemble -
    sensitivity_ensemble = Array{Array{Float64,2},1}()

    # main loop -
    for time_step_index = 1:number_of_time_steps

        # initialize sensitivity array for this time point -
        sensitivity_array = zeros(number_of_states,number_of_parameters)

        # get the time_index (row in the state array)
        simulation_time_index = time_index_array[time_step_index]

        # process each parameter -
        for parameter_index = 1:number_of_parameters

            # get parameter path -
            parameter_path = reverse(splitpath(parameter_path_array[parameter_index]))

            # what is the nominal value of this parameter?
            nominal_parameter_value = parse(Float64, get_parameter_value(parameter_path, model_dictionary))

            # load the data for this parameter -
            # nominal -
            nominal_file_path = joinpath(simulation_data_path,"simulation_nominal.dat")
            nominal_data_array = readdlm(nominal_file_path)

            # up file -
            up_file_path = joinpath(simulation_data_path, "simulation_up_P$(parameter_index).dat")
            up_data_array = readdlm(up_file_path)

            # down file -
            down_file_path = joinpath(simulation_data_path, "simulation_down_P$(parameter_index).dat")
            down_data_array = readdlm(down_file_path)

            for state_index = 1:number_of_states

                # get the state index -
                simulation_state_index = state_index_array[state_index]

                # what is the nominal state?
                nominal_state_value = nominal_data_array[simulation_time_index, (simulation_state_index+1)]

                # at this time point, calculate the difference in state -
                delta_state = up_data_array[simulation_time_index, (simulation_state_index+1)] - down_data_array[simulation_time_index, (simulation_state_index+1)]
                delta_parameter = 2*(0.05)*nominal_parameter_value

                # compute the sensitivity coeffient -
                scale_factor = (nominal_parameter_value/nominal_state_value)
                sensitivity_coefficient = scale_factor*(delta_state/delta_parameter)

                println("delta_state = $(delta_state)")

                # package in the array -
                sensitivity_array[state_index,parameter_index] = sensitivity_coefficient
            end
        end

        # push -
        push!(sensitivity_ensemble,sensitivity_array)

        # dump -
        sensitivity_output_file_path = joinpath(sensitivity_output_dir, "sensitivity_T$(time_step_index).dat")
        writedlm(sensitivity_output_file_path, sensitivity_array)
    end
end

function time_average_senstivity_array(time_index_array, simulation_data_path, sensitivity_output_dir)

    # load the nominal data  -
    nominal_file_path = joinpath(simulation_data_path,"simulation_nominal.dat")
    nominal_data_array = readdlm(nominal_file_path)

    # load sensitivity file -
    sensitivity_array_path = joinpath(sensitivity_output_dir, "sensitivity_T1.dat")
    sensitivity_array = readdlm(sensitivity_array_path)
    (number_of_states, number_of_parameters) = size(sensitivity_array)

    # how many time steps are we going to simulate over?
    number_of_time_steps = length(time_index_array)

    # get the time array -
    simulation_time_array = nominal_data_array[:,1]

    # initialize -
    time_integrated_sensitivity_array = zeros(number_of_states, number_of_parameters)

    # loop -
    for parameter_index = 1:number_of_parameters
        for state_index = 1:number_of_states

            # initialize -
            local_state_array = Float64[]
            local_time_array = Float64[]
            for time_index = 1:number_of_time_steps

                # load -
                sensitivity_array_path = joinpath(sensitivity_output_dir, "sensitivity_T$(time_index).dat")
                sensitivity_array = readdlm(sensitivity_array_path)

                # grab -
                push!(local_time_array, simulation_time_array[time_index])
                push!(local_state_array, sensitivity_array[state_index, parameter_index])
            end

            # integrate -
            delta_time = simulation_time_array[time_index_array[end]] - simulation_time_array[time_index_array[1]]
            integrated_value = (1/delta_time)*trapz(local_time_array, abs.(local_state_array))

            # cache -
            time_integrated_sensitivity_array[state_index, parameter_index] = integrated_value
        end
    end

    # return -
    return time_integrated_sensitivity_array
end

function trapz(x::Array{Float64}, y::Array{Float64})

    # Trapezoidal integration rule
    local n = length(x)
    if (length(y) != n)
        error("Vectors 'x', 'y' must be of same length")
    end

    r = 0
    if n == 1; return r; end
    for i in 2:n
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    #= correction -h^2/12 * (f'(b) - f'(a))
    ha = x[2] - x[1]
    he = x[end] - x[end-1]
    ra = (y[2] - y[1]) / ha
    re = (y[end] - y[end-1]) / he
    r/2 - ha*he/12 * (re - ra)
    =#
    return r/2
end
