function calculate_discrete_steady_state(ic_array, data_dictionary)

    # grab the step size -
    time_step_size = data_dictionary[:time_step_size]
    time_start = 0.0
    time_stop = 100.0

    # main loop -
    is_ok_to_stop = false
    steady_state_vector = zeros(1)
    max_number_of_iterations = 100
    counter = 1
    while (is_ok_to_stop == false)

        # run the model for a bit -
        (T1,soln_array_1) = GRNDiscreteDynamicSolve((time_start, time_stop, time_step_size), data_dictionary)

        # run for one more step -
        data_dictionary[:initial_condition_array] = soln_array_1[end,:]
        (T2,soln_array_2) = GRNDiscreteDynamicSolve((T1[end], T1[end]+1.0, time_step_size), data_dictionary)

        # compute the differnce -
        residual = soln_array_1[end,:] - soln_array_2[end,:]

        # error -
        error = transpose(residual)*residual

        # check error -
        if error < 0.001 || (counter > max_number_of_iterations)
            is_ok_to_stop = true
            steady_state_vector = soln_array_2[end,:]
        else
            time_start = T2[end]
            time_stop = T2[end] + 100.0
        end

        counter = counter + 1
    end

    return steady_state_vector
end

function GRNSteadyStateSolve(data_dictionary::Dict{Symbol,Any})

    # grab the initial_condition_array from the data_dictionary -
    initial_condition_array = data_dictionary[:initial_condition_array]

    # depending upon the type of problem, we will use a different steady-state routine.
    problem_type_flag = data_dictionary[:problem_type_flag]
    if problem_type_flag == :general_dynamic

        # calculate the steady-state soln -
        steady_state_prob = SteadyStateProblem(balances, initial_condition_array, data_dictionary)
        steady_state_soln  = solve(steady_state_prob, SSRootfind())

        # return -
        return steady_state_soln.u

    elseif problem_type_flag == :discrete_dynamic
        return calculate_discrete_steady_state(initial_condition_array, data_dictionary)
    else
        # Oooops - throw an error, no problem type flag
        throw(UndefVarError(:problem_type_flag))
    end
end

function GRNDiscreteDynamicSolve(time_span::Tuple{Float64, Float64, Float64}, data_dictionary::Dict{Symbol,Any})

    # Check -
    # TODO: do we have a legit time_span and data_dictionary?
    # ...
    # -------------------------------------------------------------------------- #
    # setup the problem -

    # setup the simulation time vector -
    simulation_time_array = collect(range(time_span[1], time_span[2], step=time_span[3]))

    # grab the initial_condition_array from the data_dictionary -
    initial_condition_array = data_dictionary[:initial_condition_array]

    # what is my system size?
    number_of_states = data_dictionary[:number_of_states]
    number_of_time_steps = length(simulation_time_array)

    # initalize -
    X = zeros(number_of_time_steps, number_of_states)

    # put the IC into the X array -
    for index = 1:number_of_states
        X[1,index] = initial_condition_array[index]
    end

    # main simulation loop -
    for time_value_index = 1:number_of_time_steps - 1

        # get current X -
        current_state_array = X[time_value_index,:]

        # get current time -
        current_time_value = simulation_time_array[time_value_index]

        # calculate the next state -
        next_state_array = discrete_dynamic_balances(current_time_value, current_state_array, data_dictionary)

        # package -
        for index = 1:number_of_states
            X[time_value_index+1,index] = next_state_array[index]
        end
    end

    # return -
    return (simulation_time_array, X)
end

function GRNDynamicSolve(time_span::Tuple{Float64,Float64}, data_dictionary::Dict{Symbol,Any})

    # Check -
    # TODO: do we have a legit time_span and data_dictionary?
    # ...

    # -------------------------------------------------------------------------- #
    # setup the problem -
    # grab the initial_condition_array from the data_dictionary -
    initial_condition_array = data_dictionary[:initial_condition_array]

    # build problem object -
    problem_object = ODEProblem(balances, initial_condition_array, time_span, data_dictionary)

    # solve -
    solution = solve(problem_object)

    # pull solution apart -
    T = solution.t

    # initialize the state array -
    number_of_times_steps = length(T)
    number_of_states = length(initial_condition_array)
    X = zeros(number_of_times_steps,number_of_states)
    for step_index=1:number_of_times_steps

        # grab the solution 0
        soln_array = solution.u[step_index]
        for state_index = 1:number_of_states
            X[step_index, state_index] = soln_array[state_index]
        end
    end

    # return -
    return (T,X)
    # -------------------------------------------------------------------------- #
end
