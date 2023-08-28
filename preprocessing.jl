include("/Users/Meilame/Desktop/MIT/utils.jl")

"""
    OnD_origins(tsnetwork, max_waiting, departure_times)
    Given a threshold on the maximum waiting time, returns Dict p -> Set of possible origin nodes
"""
function OnD_origins(tsnetwork, max_waiting, departure_times)
    result = Dict()
    tstep = tsnetwork.tstep
    for p in 1:tsnetwork.numcust
        sheesh = [tsnetwork.nodeid[p, t] for t in departure_times[p]:tstep:min(max_waiting, tsnetwork.horizon)] 
        result[p] = sheesh
    end

    return result
end

"""
    OnD_destinations(tsnetwork, max_arrival_times, departure_times)
    Given a threshold on the maximum arrival time, returns Dict p -> Set of possible destination nodes
"""
function OnD_destinations(tsnetwork, max_arrival_times, departure_times)
    result = Dict()
    for p in 1:tsnetwork.numcust
        t_to_dest = determine_distance(tsnetwork, p, p+tsnetwork.numcust)
        sheesh = [tsnetwork.nodeid[p + tsnetwork.numcust, t] for t in departure_times[p] + t_to_dest:min(tsnetwork.horizon, max_arrival_times[p])]
        result[p] = sheesh
    end

    return result
end


"""
    buildRoutes(transit_arc_list)
    Build ALL the possible routes (Number = number of transit arcs * 5)
    Returns a dictionary routeid[arc, type] --> route_id, and routedesc[route_id] --> (arc, type), and the number of routes
"""
function buildRoutes(transit_arc_list)
    routeid = Dict()
    routedesc = Dict()
    index = 1
    println("Building routes...")
    for arc in transit_arc_list
        for type in 1:5
            if type != 2
                routeid[arc, type] = index
                routedesc[index] = (arc, type)
                index += 1
            end
        end
    end

    #Type 2 has no transit arc (empty). Arc number 0 and route id 0 will materialize type 2 route
    routeid[0, 2] = 0
    routedesc[0] = (0, 2)
    if index == length(collect(keys(routedesc)))
        println("Done. Total number of routes = ", index)
    end
    return routeid, routedesc, index
end


"""
    determineRouteWalkingDistance(tsnetwork, numroutes, routedesc, times_to_station, dist_alight_dest)
    Returns Dict[p, route_id] => Walking distance for passenger p if assigned to route route_id
"""
function determineRouteWalkingDistance(tsnetwork, numroutes, routedesc, times_to_station, dist_alight_dest)
    res = Dict()
    for p in 1:tsnetwork.numcust
        for route_id in 1:numroutes-1
            arc, type = routedesc[route_id]
            boarding_loc = tsnetwork.nodedesc[tsnetwork.arcdesc[arc][1]][1]
            alighting_loc = tsnetwork.nodedesc[tsnetwork.arcdesc[arc][2]][1]
            if type == 1
                res[p, route_id] = times_to_station[p, boarding_loc] + dist_alight_dest[p, alighting_loc]
            end
            if type == 3
                res[p, route_id] = dist_alight_dest[p, alighting_loc]
            end
            if type == 4
                res[p, route_id] = times_to_station[p, boarding_loc]
            end
            if type == 5
                res[p, route_id] = 0
            end
        end
        res[p, 0] = 0
    end
    return res
end

function determineBoardingTime(tsnetwork, arc_id)
    n1, n2, dist = tsnetwork.arcdesc[arc_id]
    return tsnetwork.nodedesc[n1][2]
end

#Max arrival time defined as the time we would arrive if we had walked
function maxArrivalTime(tsnetwork, walking_speed)
    dist_to_dest = distance_to_dest(tsnetwork)
    res = Dict()
    for p in 1:tsnetwork.numcust
        res[p] = dist_to_dest[p] / walking_speed
    end
    return res   
end


"""
    boardingNodeRoute(tsnetwork, numroutes, routedesc)
    Returns Dict[route_id] --> boarding node (in the OnD tsnetwork)
"""
function boardingNodeRoute(tsnetwork, numroutes, routedesc)
    res = Dict()
    for route in 1:numroutes-1
        transit_arc, type = routedesc[route]
        boarding_node, alight_node, transit_distance = tsnetwork.arcdesc[transit_arc]
        res[route] = boarding_node
    end
    return res
end

#Returns Dict[node] --> Set of routes which boarding node is node
function reverseBoardingNodeRoute(boarding_node_route)
    res = Dict()
    for b_node in collect(values(boarding_node_route))
        res[b_node] = Set()
    end
    for route in keys(boarding_node_route)
        push!(res[boarding_node_route[route]], route)
    end
    return res
end

#Returns Dict[route] --> Alight node (index in the OnD tsnetwork)
function alightNodeRoute(tsnetwork, numroutes, routedesc)
    res = Dict()
    for route in 1:numroutes-1
        transit_arc, type = routedesc[route]
        boarding_node, alight_node, transit_distance = tsnetwork.arcdesc[transit_arc]
        res[route] = alight_node
    end
    return res
end

#Returns Dict[node] --> Set of routes which alight node is node
function reverseAlightNodeRoute(alight_node_route)
    res = Dict()
    for a_node in collect(values(alight_node_route))
        res[a_node] = Set()
    end
    for route in keys(alight_node_route)
        push!(res[alight_node_route[route]], route)
    end
    return res
end

#Returns Dict[route] --> Boarding station ∈ {2*numcust + 1, ..., 2*numcust + numstations}
function boardingStationRoute(tsnetwork, numroutes, routedesc)
    temp = boardingNodeRoute(tsnetwork, numroutes, routedesc)
    res = Dict()
    for r in 1:numroutes-1
        node = temp[r]
        res[r] = tsnetwork.nodedesc[node][1]
    end
    return res
end

#Returns Dict[route] --> Alighting station ∈ {2*numcust + 1, ..., 2*numcust + numstations}
function alightingStationRoute(tsnetwork, numroutes, routedesc)
    res = Dict()
    temp = alightNodeRoute(tsnetwork, numroutes, routedesc)
    for r in 1:numroutes-1
        node = temp[r]
        res[r] = tsnetwork.nodedesc[node][1]
    end
    return res
end


"""
    arrivalTimesType13(tsnetwork, numroutes, numstations, route, routedesc, walking_speed, horizon, p)
    For a given p, returns for routes of type 1 and 3 the implied arrival time at the destination (basically, we walk from alighting station to destination)
    Returns Int
"""
function arrivalTimesType13(tsnetwork, numroutes, numstations, route, routedesc, walking_speed, horizon, p)
    arc_times = arc_arrival_times(tsnetwork)
    dist_to_dest = distance_alight_to_dest(tsnetwork, numstations)
    alight_station = alightingStationRoute(tsnetwork, numroutes, routedesc)[route]

    transit_arc = routedesc[route][1]
    transit_arrival = arc_times[transit_arc]
    walking_time = dist_to_dest[p, alight_station] / walking_speed
    return min(horizon, transit_arrival + walking_time)
end


"""
    route_boarding_alighting_time(tsnetwork, boarding_node_route, alight_node_route, numroutes)
    boarding_node_route : Dict[route] --> boarding node
    alight_node_route : Dict[route] --> alight node
    Returns Dict[route] --> boarding time, Dict[route] --> alighting time
"""
function route_boarding_alighting_time(tsnetwork, boarding_node_route, alight_node_route, numroutes)
    board_time = Dict()
    alight_time = Dict()
    for r in 1:numroutes-1
        b_node = boarding_node_route[r]
        a_node = alight_node_route[r]
        board_time[r] = tsnetwork.nodedesc[b_node][2]
        alight_time[r] = tsnetwork.nodedesc[a_node][2]
    end
    return board_time, alight_time
end



"""
    preprocessing(ond_tsnetwork, transit_tsnetwork, num_stations, max_walking, max_waiting, walking_speed, departure_times)
    Preprocessing algorithm
"""
function preprocessing(ond_tsnetwork, transit_tsnetwork, num_stations, max_walking, max_waiting, walking_speed, departure_times)
"""
BASICALLY, A ROUTE IS (TRANSIT_ARC, TYPE)

R1, R2 .. R5 such as R1[p] = [list of route_id], with a mapping Dict[route_id] --> transit_arc
and 
R such as R[p] = union(R1[p], ..., R5[p])

This function should return R, R1 .. R5, Npo, Npb, Npa, Npd adjusted
"""

    println("Beginning preprocessing...")
    transit_arc_list, boarding_nodes = transit_arclist_OnD(ond_tsnetwork, transit_tsnetwork)
    routeid, routedesc, numroutes = buildRoutes(transit_arc_list)
    times_to_station = distance_to_station(ond_tsnetwork, num_stations)
    walking_times_to_station = walking_t_to_station(ond_tsnetwork, num_stations, walking_speed)
    distance_alight_dest = distance_alight_to_dest(ond_tsnetwork, num_stations)
    max_arrival_times = maxArrivalTime(ond_tsnetwork, walking_speed) #we should arrive before the time we would arrive if we had walked ; Dict p -> max_time
    walking_dict = determineRouteWalkingDistance(ond_tsnetwork, numroutes, routedesc, times_to_station, distance_alight_dest)

    # --- Init dictionaries --- #
    R1, R2, R3, R4, R5, R = Dict(), Dict(), Dict(), Dict(), Dict(), Dict()
    R_dicts = [R1, R2, R3, R4, R5]
    Npo, Npb, Npa, Npd = Dict(), Dict(), Dict(), Dict()
    for p in 1:ond_tsnetwork.numcust
        R1[p], R2[p], R3[p], R4[p], R5[p] = [],[0],[],[],[]
        Npo[p], Npb[p], Npa[p], Npd[p] = Set(), Set(),Set(),Set()
    end
    #---------------------------#
    temp = Set() #counting purpose
    push!(temp, 0)
    for route in 1:numroutes-1
        transit_arc, type = routedesc[route]
        R_type = R_dicts[type]
        t_boarding = determineBoardingTime(ond_tsnetwork, transit_arc)
        boarding_node, alight_node, transit_distance = ond_tsnetwork.arcdesc[transit_arc]
        boarding_station = ond_tsnetwork.nodedesc[boarding_node][1]
        alighting_station = ond_tsnetwork.nodedesc[alight_node][1]

        for p in 1:ond_tsnetwork.numcust
            tto_station = times_to_station[p, boarding_station]
            walking_time_to_station = walking_times_to_station[p, boarding_station]
            if walking_dict[p, route] <= max_walking
                if type == 1 || type == 4
                    if departure_times[p] + walking_time_to_station <= t_boarding <= departure_times[p] +  walking_time_to_station + max_waiting
                            push!(R_type[p], route)
                            push!(temp, route)
                            if type == 4
                                push!(Npa[p], alight_node)
                            end
                    end
                end
                if type == 3 || type == 5       
                    if departure_times[p] + tto_station <= t_boarding <= tto_station + max_waiting
                        push!(R_type[p], route)
                        push!(temp, route)
                        if type == 3
                            push!(Npb[p], boarding_node)
                        end
                        if type == 5
                            push!(Npb[p], boarding_node)
                            push!(Npa[p], alight_node)
                        end  
                    end  
                end
            end 
        end
    end
    
    #Building R 
    for p in 1:ond_tsnetwork.numcust
        R[p] = union(R1[p], R2[p], R3[p], R4[p], R5[p]) #We note that type 2 route always feasible whatever the constraints, so route n° 0 will always be there
    end
    Npo = OnD_origins(ond_tsnetwork, max_waiting, departure_times)
    Npd = OnD_destinations(ond_tsnetwork, max_arrival_times, departure_times)

    println("PREPROCESSING DONE. Updated number of routes: ", length(temp))

    return R, R1, R2, R3, R4, R5, Npo, Npb, Npa, Npd, routeid, routedesc, numroutes, walking_dict

end



