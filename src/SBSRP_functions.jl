using Plots, LightGraphs, GraphRecipes, TravelingSalesmanExact, GLPK

function greedy_algo_SBRP(c, initial_inventory)

    tour       = [1]
    inventory  = deepcopy(initial_inventory)
    bikes_held = 0

    while !all(inventory .== 0)

        potential_cities = c[:, tour[end]]

        # decide where to go next based on current number of bikes
        if bikes_held > 0
            potential_cities[findall(inventory .>= 0)] .= Inf
        else
            potential_cities[findall(inventory .<= 0)] .= Inf
        end
        
        next_station = findmin(potential_cities)[2]

        # decide how many bikes leave and how many take
        leave = min(0, (inventory[next_station] + bikes_held))
        take  = max(0, (inventory[next_station] + bikes_held))

        inventory[next_station] = leave
        bikes_held = take

        push!(tour, next_station)
        println("station: ", next_station, " | inventory: ", inventory)
    end
    return [inventory, tour]
end

# for a given tour and initial_inventory of equal length check what final inventory looks
function check_inventory(tour, initial_inventory)

    inv  = deepcopy(initial_inventory)
    bikes_held = 0

    for vertex in tour

        leave = min(0, (inv[vertex] + bikes_held))
        take  = max(0, (inv[vertex] + bikes_held))

        inv[vertex] = leave
        bikes_held = take

    end

    return [inv, bikes_held]
end

function proba_algo(c, initial_inventory)

    # find solution as if its TSP
    set_default_optimizer!(GLPK.Optimizer)
    tsp_solved      = get_optimal_tour(c)[1]
    proposition     = deepcopy(tsp_solved)

    # check what is the balance at bike stations affter TSP solution
    invent_solution = check_inventory(tsp_solved, initial_inventory)
    inventory       = []

    # for each station with demand
    for demand_station in findall(invent_solution .!= 0)

        # get top 5 nearest stations of those with demand
        close_stations = c[:, demand_station]
        close_stations = sortperm(close_stations)[(end-5):end]

        # for each top 5 nearest stations with those with demand
        for near_station in reverse(close_stations)

            # swap nearest with the one with demand
            swap_proposition = deepcopy(proposition)
            swap_proposition[near_station] = swap_proposition[demand_station]

            # check what is the balance at stations after swap
            balance = check_inventory(swap_proposition, initial_inventory)

            # if every station is satisfied then break 
            if all(balance .== 0)
                feasible_proposition = swap_proposition
                break
            elseif balance[demand_station] == 0
                proposition = deepcopy(swap_proposition)
                break
            end
        end

        if @isdefined feasible_proposition
            break
        end
    end

    if @isdefined feasible_proposition
        return feasible_proposition
    end
end


function adjust_TSP(dist_mat, bikes)

    # find solution as if its TSP
    set_default_optimizer!(GLPK.Optimizer)
    tsp_solved      = get_optimal_tour(dist_mat)[1]
    proposition     = deepcopy(tsp_solved)

    # check what is the balance at bike stations affter TSP solution
    invent_solution = check_inventory(tsp_solved, bikes)[1]
    in_demand = findall(invent_solution .!= 0)
    inventory = deepcopy(invent_solution)

    while !all(inventory .== 0)

        # get distances from station with demand
        nearest_station = dist_mat[:, in_demand[1]]
        # eliminate those with no excess
        nearest_station[bikes .<= 0] .= Inf
        # eliminate visited
        visited = proposition[1:findall(proposition .== in_demand[1])[1]]
        nearest_station[visited] .= Inf
        # which index in proposition is the best candidate
        nearest_station = findall(proposition .== findmin(nearest_station)[2])[1]
        # which indexes of proposition to swap
        which_prop_swap = [findall(proposition .== in_demand[1])[1], nearest_station]
        # swap them
        proposition[which_prop_swap[2]], proposition[which_prop_swap[1]] = proposition[which_prop_swap[1]], proposition[which_prop_swap[2]]
        # check how has demand changed
        new_inventory = check_inventory(proposition, bikes)[1]
        
        if length(findall(invent_solution .!= 0)) > length(findall(new_inventory .!= 0))
            inventory = new_inventory
        end

        in_demand = findall(new_inventory .!= 0)

        if length(in_demand) == 0
            global solved_tour = proposition
        end
    end
    return solved_tour
end

function solve_asif_TSP(c, bikes)

    # find solution as if its TSP
    set_default_optimizer!(GLPK.Optimizer)
    tsp_solved      = get_optimal_tour(c)[1]
    proposition     = deepcopy(tsp_solved)

    # check what is the balance at bike stations affter TSP solution
    invent_solution = check_inventory(tsp_solved, bikes)[1]
    in_demand = findall(invent_solution .!= 0)

    while true
        
        nearest_station = c[tsp_solved[end], :]
        nearest_station = findmin(nearest_station[in_demand])[2]
        nearest_station = in_demand[nearest_station]

        push!(tsp_solved, nearest_station)
        splice!(in_demand, findall(in_demand .== nearest_station)[1])

        if length(in_demand) == 0
            break
        end
    end
    return tsp_solved
end
