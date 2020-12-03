using JuMP, GLPK
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
    N    = size(x)[1]
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

# euclidean distance (function with 2 methods)

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

