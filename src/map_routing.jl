using Random, XLSX, Parameters, OpenStreetMapX, LightGraphs, Colors, Plots, GraphRecipes, TravelingSalesmanExact, GLPK, OpenStreetMapXPlot

cd("C:\\Users\\HP\\Documents\\julia\\BikeSharingRebalancing.jl")
include("SBSRP_functions.jl")
Random.seed!(0)

stations = XLSX.readxlsx("bike_sharing_stations.xlsx")["Sheet1"][:]
torun = get_map_data("torun.osm", use_cache = false, trim_to_connected_graph = true)

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


# rozwiązywanie
tour_greedy = greedy_algo_SBRP(torun_BSS_dist_mat, torun_bikes)[2]
tour_asif = solve_asif_TSP(torun_BSS_dist_mat, torun_bikes)
tour_adjust = adjust_TSP(torun_BSS_dist_mat, torun_bikes)

# ile tras wyznaczono
tour_greedy |> get_graph_matrix |> sum
tour_asif |> get_graph_matrix |> sum
tour_adjust |> get_graph_matrix |> sum

# czy wszystkie uzupelnily stacje?
all(check_inventory(tour_adjust, torun_bikes)[1] .== 0)
all(check_inventory(tour_adjust, torun_bikes)[1] .== 0)
all(check_inventory(tour_adjust, torun_bikes)[1] .== 0)

# ile wynosi funkcja celu (długość trasy)
sum(get_graph_matrix(tour_greedy) .* torun_BSS_dist_mat)
sum(get_graph_matrix(tour_asif) .* torun_BSS_dist_mat)
sum(get_graph_matrix(tour_adjust) .* torun_BSS_dist_mat)

plot(stations[:,1],stations[:,2], seriestype = :scatter, markersize = 3, color = :black, ticks = false,
     axis = false, legend=:none, label = "miasta", series_annotation = text.(torun_bikes, :bottom))
     plot!(stations[:,1][tour_adjust, 1], stations[:,2][tour_adjust, 1], linecolor = :red, label = "ścieżka")



Random.seed!(0);
pointA = station_nodes[10]
pointB = station_nodes[39]
sr = shortest_route(torun, pointA, pointB)[1]
Plots.gr()
p = OpenStreetMapXPlot.plotmap(torun,width=1200,height=800);
addroute!(p, torun, sr;route_color="red", fontsize = 10);

tour = Int64[]
for station in 1:(length(station_nodes)-1)
    global tour = vcat(tour, shortest_route(torun, station_nodes[station], station_nodes[station + 1])[1])
end


tour = [push!(tour, shortest_route(torun, station_nodes[a], station_nodes[a + 1])[1]) for a = 1:(length(station_nodes)-1)]

p = OpenStreetMapXPlot.plotmap(torun,width=1200,height=800)
addroute!(p, torun, tour;route_color="red", fontsize = 10)

for point in 1:length(station_nodes) 
    addroute!(p, torun, )
end

shortest_route(torun, station_nodes[1], station_nodes[11])[1][:]

s = []
push!(s, [2,3,4,2])
sr
Int64{tour}
