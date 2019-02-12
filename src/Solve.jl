function GRNSteadyStateSolve(data_dictionary::Dict{Symbol,Any})

    # grab the initial_condition_array from the data_dictionary -
    initial_condition_array = data_dictionary[:initial_condition_array]

    # calculate the steady-state soln -
    steady_state_prob = SteadyStateProblem(balances, initial_condition_array, data_dictionary)
    steady_state_soln  = solve(steady_state_prob, SSRootfind())

    # return -
    return steady_state_soln.u
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
