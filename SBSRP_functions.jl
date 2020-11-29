using Plots, LightGraphs, GraphRecipes, TravelingSalesmanExact, GLPK

# constraints

# each vertice can only have two edges
function double_edge(x)
    N = size(x)[1]
    [sum(x[city,:]) == 1 for city=1:N]
    [sum(x[:,city]) == 1 for city=1:N]
end

# no city shall be connected to itself
function self_connected(x)
    N = size(x)[1]
    [x[city,city] == 0 for city=1:N]
end

# if there is an edge f to t then there can be no t to f edge
function come_back(x)
    N = size(x)[1]
    [x[f,t] + x[t,f] <= 1 for t = 1:N, f = 1:N]
end

# checks if there is only one tour i.e. no subtours
function single_tour(x)
    N = size(x)[1]
    tour = int[1]
    while true
        city_loc = findmax(x[tour[end], 1:N])
        if city_loc == tour[1]
            break
        else
            push!(tour, city_loc)
        end
    end

    if length(tour) <= N
        false
    end
    true
end

# imput a vector of cities in order of tour and get matrix of tour 
function get_graph_matrix(tour)
    N = maximum(tour)
    g = zeros(N,N)
    y_cord = tour
    x_cord = tour[1:end-1]
    pushfirst!(x_cord, tour[end])
    
    for cord in 1:length(x_cord)
        g[x_cord[cord], y_cord[cord]] += 1
    end
    g
end

# euclidean distance
function euc_distance(p_1, p_2, q_1, q_2)
    sqrt.((p_1 .- q_1).^2 .+ (p_2 .- q_2).^2)
end

# make distance_matrix. Another method for different arguments.
function euc_distance(x_coord, y_coord)
    dist_mat = zeros(length(x_coord),length(x_coord))
    for x in 1:N 
        for y in 1:N 
            dist_mat[x,y] = sqrt.((x_coord[x] .- x_coord[y]).^2 .+ (y_coord[x] .- y_coord[y]).^2)
        end
    end
    return dist_mat
end

# algorithm solving by choosing nearest city until the last one
function greedy_algo(c)
    tour = [1]
    while true
        potential_cities = c[:, tour[end]]
        potential_cities[tour] .= Inf
        city = findmin(potential_cities)[2]
        if city == tour[1]
            break
        else
            push!(tour, city)
        end
        println("Next city: ", city)
    end
    println("found path")
    splice!(tour, 1)
    push!(tour, 1)
    tour
end

# simulate inbalanced bike sharing stations
# generate random vector of n_city length that sums to 0 and have 
# max equal to max_inventory and min equal to -max_inventory
function inbalance(n_city, max_inventory)

    inventory = rand(-max_inventory:max_inventory, n_city)

    while sum(inventory) != 0
        which_change = rand(1:length(inventory))
        add =  sum(inventory) > 0 ? -1 : 1

        if inventory[which_change] == max_inventory
            inventory[which_change] = inventory[which_change] - 1
        elseif inventory[which_change] == -max_inventory
            inventory[which_change] = inventory[which_change] + 1
        else
            inventory[which_change] = inventory[which_change] + add
        end
    end
    return inventory
end


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

function adjust_TSP(c, initial_inventory)

    # find solution as if its TSP
    set_default_optimizer!(GLPK.Optimizer)
    tsp_solved      = get_optimal_tour(c)[1]
    proposition     = deepcopy(tsp_solved)

    # check what is the balance at bike stations affter TSP solution
    invent_solution = check_inventory(tsp_solved, initial_inventory)[1]
    in_demand = findall(invent_solution .!= 0)
    inventory = []

    while !all(inventory .== 0)
        nearest_stations = c[:, in_demand[1]]
        nearest_stations = nearest_stations[findall(initial_inventory .> 0)]
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

        nearest_station = dist_mat[:, in_demand[1]]
        nearest_station[findall(bikes .<= 0)] .= Inf
        visited = proposition[1:findall(proposition .== in_demand[1])[1]]
        nearest_station[visited] .= Inf
        nearest_station = findall(proposition .== findmin(nearest_station)[2])[1]

        proposition[in_demand[1]], proposition[nearest_station] = proposition[nearest_station], proposition[in_demand[1]]

        new_inventory = check_inventory(proposition, bikes)[1]
        
        if length(findall(invent_solution .!= 0)) > length(findall(new_inventory .!= 0))
            inventory = new_inventory
        end

        in_demand = findall(new_inventory .!= 0)
    end
    return [inventory, proposition]
end

function adjust_TSP2(dist_mat, bikes)

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
