using XLSX, Random, OpenStreetMapX, Colors, Plots, LightGraphs, GraphRecipes, TravelingSalesmanExact, GLPK, OpenStreetMapXPlot

include("src\\SBSRP_functions.jl")
include("src\\misc_functions.jl")
Random.seed!(0)

stations = XLSX.readxlsx("geo\\bike_sharing_stations.xlsx")["Sheet1"][:]
torun = get_map_data("geo\\torun.osm", use_cache = false, trim_to_connected_graph = true)

n_stations = size(stations)[1]
torun_BSS_dist_mat = zeros(n_stations, n_stations)
station_nodes = [point_to_nodes((stations[s,1], stations[s,2]), torun) for s = 1:size(stations)[1]] 

for station_a in 1:n_stations
    for station_b in 1:n_stations 

        if station_b == station_a 
            torun_BSS_dist_mat[station_a, station_b] = 0
        elseif torun_BSS_dist_mat[station_b, station_a] != 0
            torun_BSS_dist_mat[station_a, station_b] = torun_BSS_dist_mat[station_b, station_a]
        else
            torun_BSS_dist_mat[station_a, station_b] = shortest_route(torun, station_nodes[station_a], station_nodes[station_b])[2]
        end   

        println("done: ", round((station_a * (n_stations-1) + station_b)/(n_stations^2) * 100, digits = 3), " %")      

    end
end

torun_bikes = inbalance(size(torun_BSS_dist_mat)[1], 7)


# solving 
tour_greedy = greedy_algo_SBRP(torun_BSS_dist_mat, torun_bikes)[2]
tour_asif = relax_greedy_SBSRP(torun_BSS_dist_mat, torun_bikes)
tour_adjust = permute_relax_SBSRP(torun_BSS_dist_mat, torun_bikes)

# through how many stations proposed tour goes?
tour_greedy |> get_graph_matrix |> sum
tour_asif |> get_graph_matrix |> sum
tour_adjust |> get_graph_matrix |> sum

# are stations full?
all(check_inventory(tour_adjust, torun_bikes)[1] .== 0)
all(check_inventory(tour_adjust, torun_bikes)[1] .== 0)
all(check_inventory(tour_adjust, torun_bikes)[1] .== 0)

# Whats the value of objective function (length of tour)
sum(get_graph_matrix(tour_greedy) .* torun_BSS_dist_mat[1:maximum(tour_greedy), 1:maximum(tour_greedy)]) # greedy algo omit stations with no demand so tour graph is smaller
sum(get_graph_matrix(tour_asif) .* torun_BSS_dist_mat)
sum(get_graph_matrix(tour_adjust) .* torun_BSS_dist_mat)

plot(stations[:,1],stations[:,2], seriestype = :scatter, markersize = 3, color = :black, ticks = false,
     axis = false, series_annotation = text.(torun_bikes, :bottom), label = "stacje")


tour_animate = deepcopy(tour_adjust)

BSRP_animation = @animate for station in 1:length(tour_animate)

    tour_loop = tour_animate[1:station]
    station_d = check_inventory(tour_loop, torun_bikes)[1]

    plot(stations[:,1],stations[:,2], seriestype = :scatter, markersize = 3,
    color = :black, ticks = false, 
    series_annotation = text.(station_d, :bottom), label = "Nodes/stations")
    plot!(stations[tour_loop,1], stations[tour_loop, 2], label = "Vertices/tour")

end

gif(BSRP_animation, fps = 2,  "doc\\media\\adjusting_algo.gif")





