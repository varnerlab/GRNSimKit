function calculate_discrete_steady_state(ic_array, data_dictionary)

    # setup the optimization problem -
    objective_function(x) = discrete_steady_state_balance_residual(x,data_dictionary)

    # get the lower, and upper_bounds for the states -
    lower_bound_array = data_dictionary[:lower_bound_array]
    upper_bound_array = data_dictionary[:upper_bound_array]

    # make a call to the optim package -
    result = optimize(objective_function,lower_bound_array, upper_bound_array, ic_array,Fminbox(LBFGS()))

    # get the extent -
    steady_state_vector = Optim.minimizer(result)

    # return -
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
