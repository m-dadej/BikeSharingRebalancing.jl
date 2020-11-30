
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
