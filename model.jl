include("utils.jl")
include("preprocessing.jl")
using JuMP, Gurobi

"""
The Model
"""
function model(ond_tsnetwork, numcust, numstations, numzones, horizon, K, W, R, R1, R2, R3, R4, R5, Npo, Npb, Npa, Npd, walking_speed, routedesc, numroutes, TIME_COSTS, MOVING_COST, WALKING_COST)
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", 1)
    time_limit = 1800
    set_time_limit_sec(model, time_limit)

    #---------------------UTILS----------------------------#
preprocessing_t = @elapsed begin
    #Zones separation
    I = define_zones(ond_tsnetwork)

    other_nodes = Dict()
    #to be optimized...
    for p in 1:numcust
        other_nodes[p] = []
        for node in 1:ond_tsnetwork.numnodes
            if !(node in union(values(Npo[p]), values(Npb[p]), values(Npa[p]), values(Npd[p])))
                push!(other_nodes[p], node)
            end
        end
    end


    walking_times_to_station = walking_t_to_station(ond_tsnetwork, numstations, walking_speed) #Dict[p, s] --> walking time from origin of p to station s

    moving_arcs, waiting_arcs = partition_arcs(ond_tsnetwork) 
    flow_balance_nodes = flowBalanceNodes(ond_tsnetwork) #for vehicles

    dist_alight_to_dest = distance_alight_to_dest(ond_tsnetwork, numstations) #Dict[p, s] --> distance from station s to destination of p
    station_waiting_arcs = stationWaitingArcs(ond_tsnetwork, waiting_arcs,numstations) #arcs that are waiting arcs and that are at a station in the on-demand TSN

    f_arc_keys = union(moving_arcs, station_waiting_arcs) #basically, \mathcal{A}^{tr} 

    nontransfer_arcs = 1:ond_tsnetwork.numarcs

    routeWalkingDict = determineRouteWalkingDistance(ond_tsnetwork, numroutes, routedesc, walking_times_to_station, dist_alight_to_dest) #Dict[p, r] --> walking distance for passenger p on route r

    boarding_route = boardingNodeRoute(ond_tsnetwork, numroutes, routedesc) #Dict[route] --> boarding node
    alighting_route = alightNodeRoute(ond_tsnetwork, numroutes, routedesc) #Dict[route] --> alighting node
    reverse_boarding_route = reverseBoardingNodeRoute(boarding_route) #Dict[node] --> Set () routes that boards at n
    reverse_alighting_route = reverseAlightNodeRoute(alighting_route) #Dict[node] --> Set () routes that boards at n

    #--------------------- VARIABLES -------------------------------# 

    @variable(model, 0<= y[a in nontransfer_arcs], Int)
    @variable(model, 0 <= f[p in 1:numcust, a in f_arc_keys] <= 1)
    @variable(model, x[p in 1:numcust, r in R[p]], Bin)
    @variable(model, zO[p in 1:numcust, n in Npo[p]], Bin)
    @variable(model, zD[p in 1:numcust, n in Npd[p]], Bin)


    #--------------------- CONSTRAINTS -------------------------------#

    @constraint(model, [p in 1:numcust], sum(x[p, r] for r in R[p]) == 1)
    @constraint(model, [p in 1:numcust], sum(zO[p, n] for n in Npo[p]) == sum(x[p,r] for r in R2[p]) + sum(x[p,r] for r in R3[p]) + sum(x[p,r] for r in R5[p]))
    @constraint(model, [p in 1:numcust], sum(zD[p, n] for n in Npd[p]) == sum(x[p,r] for r in R2[p]) + sum(x[p,r] for r in R4[p]) + sum(x[p,r] for r in R5[p]))


    @constraint(model, [p in 1:numcust, n in Npo[p]], sum(f[p,a] for a in intersect(ond_tsnetwork.A_plus[n], f_arc_keys)) - sum(f[p,a] for a in intersect(ond_tsnetwork.A_minus[n], f_arc_keys)) == zO[p,n])
    @constraint(model, [p in 1:numcust, n in Npb[p]], sum(f[p,a] for a in intersect(ond_tsnetwork.A_plus[n], f_arc_keys)) - sum(f[p,a] for a in intersect(ond_tsnetwork.A_minus[n], f_arc_keys)) == - sum(x[p,r] for r in intersect(union(R3[p], R5[p]), reverse_boarding_route[n])))
    @constraint(model, [p in 1:numcust, n in Npa[p]], sum(f[p,a] for a in intersect(ond_tsnetwork.A_plus[n], f_arc_keys)) - sum(f[p,a] for a in intersect(ond_tsnetwork.A_minus[n], f_arc_keys)) == sum(x[p,r] for r in intersect(union(R4[p], R5[p]), reverse_alighting_route[n])))
    @constraint(model, [p in 1:numcust, n in Npd[p]], sum(f[p,a] for a in intersect(ond_tsnetwork.A_plus[n], f_arc_keys)) - sum(f[p,a] for a in intersect(ond_tsnetwork.A_minus[n], f_arc_keys)) == - zD[p,n])
    @constraint(model, [p in 1:numcust, n in other_nodes[p]], sum(f[p,a] for a in intersect(ond_tsnetwork.A_plus[n], f_arc_keys)) - sum(f[p,a] for a in intersect(ond_tsnetwork.A_minus[n], f_arc_keys)) ==0)



    #Vehicle constraints
    @constraint(model, [a in f_arc_keys], sum(f[p,a] for p in 1:numcust) <= K * y[a])
    @constraint(model, [n in flow_balance_nodes], sum(y[a] for a in intersect(ond_tsnetwork.A_plus[n], nontransfer_arcs)) - sum(y[a] for a in intersect(ond_tsnetwork.A_minus[n], nontransfer_arcs)) == 0)


    for i in 1:length(I)
        @constraint(model, [1:i], sum(sum(y[a] for a in intersect(ond_tsnetwork.A_plus[n], nontransfer_arcs)) for n in I[i]) == W[i])
    end


    @objective(model, Min,
    sum(TIME_COSTS[ond_tsnetwork.nodedesc[node][2]]*zO[p, node] for p in 1:numcust, node in Npo[p]) 
    + sum(TIME_COSTS[ond_tsnetwork.nodedesc[node][2]]*zD[p, node] for p in 1:numcust, node in Npd[p])
    + sum(TIME_COSTS[arrivalTimesType13(ond_tsnetwork, numroutes, numstations, r, routedesc, walking_speed, horizon, p)] * x[p,r] for p in 1:numcust, r in union(R1[p], R3[p]))
    + WALKING_COST * sum( routeWalkingDict[p, r] * x[p,r] for p in 1:numcust, r in R[p])
    + MOVING_COST * (sum(f[p,a] * ond_tsnetwork.arccost[a] for p in 1:numcust, a in moving_arcs) + sum(y[a] * ond_tsnetwork.arccost[a] for a in intersect(moving_arcs, nontransfer_arcs)))
    )  
    println("")
end
    println("OPTIMIZING MODEL...")
optimizing_t = @elapsed begin   
    optimize!(model)
    println("DONE")
end
    println("")
    println("OBJECTIVE VALUE = ", objective_value(model))

    return value.(f), value.(y), value.(x), value.(zO), value.(zD), preprocessing_t, optimizing_t
end

#--------TO VISUALIZE SOLUTIONS ONCE OPTIMIZED-------------#

"""
    f_to_arctype(f, types, tsnetwork, numstations)
    Input: optimized f, optimized x, types (Dict[p] => type), on-demand TSN, numstations
    Returns Dict[a => type] where a is an arc that is used by a passenger and type is the type assigned to this passenger
"""
function f_to_arctype(f, types, tsnetwork, numstations)
    moving_arcs, waiting_arcs = partition_arcs(tsnetwork)
    station_waiting_arcs = stationWaitingArcs(tsnetwork, waiting_arcs,numstations)
    f_arc_keys = union(moving_arcs, station_waiting_arcs)
    arclist = Dict()
    for p in 1:tsnetwork.numcust
        typ = types[p]
        for a in f_arc_keys
            if f[p,a] > 0
                arclist[a] = typ
            end
        end
    end
    return arclist
end

"""
    y_to_arclist(y, tsnetwork)
    Input: optimized y, on-demand TSN
    Returns Set of arcs that are used by a vehicle
"""
function y_to_arclist(v, tsnetwork)
    arclist = Set()
    for a in 1:tsnetwork.numarcs
        if v[a] > 0
            push!(arclist, a)
        end
    end
    return arclist
end



"""
    path(p, f, types, ond_tsnetwork, numstations)
    For a given passenger p, returns Dict[a] => type where a is an arc that is used by p and type is the type assigned to p
    It enables to draw its path on the on-demand TSN
"""
function path(p, f, types, ond_tsnetwork, numstations)
    moving_arcs, waiting_arcs = partition_arcs(ond_tsnetwork)
    station_waiting_arcs = stationWaitingArcs(ond_tsnetwork, waiting_arcs, numstations)
    f_arc_keys = union(moving_arcs, station_waiting_arcs)
    arc_to_type = Dict()
    type = types[p]
    for a in f_arc_keys
        if f[p,a] > 0
            arc_to_type[a] = type
        end
    end
    return arc_to_type
end


"""
    selectedNode(tsnetwork, z, N)
    Input: on-demand TSN, optimized z, corresponding N (O, A, B, D)
    Output: Dict[p => n] where n is the node selected for passenger p
"""
function selectedNode(tsnetwork, z, N)
    res = Dict()
    for p in 1:tsnetwork.numcust
        for n in N[p]
            if z[p,n] > 0.5
                res[p] = n
            end
        end
    end
    return res
end

#--------ANALYSIS-------------#

"""
    giveRoute(tsnetwork, x, R)
    Returns Dict[p] => r where r is the route chosen for passenger p
"""
function giveRoute(tsnetwork, x, R)
    res = Dict()
    for p in 1:tsnetwork.numcust
        for r in R[p]
            if x[p, r] > 0.5
                res[p] = r
            end
        end
    end
    return res
end


"""
    give_type_and_walk(tsnetwork, x, routedesc, R, walking_dict)
    Returns: 
    - Dict["Passenger p" (string)] => type where type is the type assigned to passenger p
    - Dict[p] => type where type is the type assigned to passenger p
    - Dict["Passenger p" (string)] => walking distance for passenger p
"""
function give_type_and_walk(tsnetwork, x, routedesc, R, walking_dict)
    chosen_routes = giveRoute(tsnetwork, x, R)
    res = Dict()
    res1 = Dict()
    walking = Dict()
    for p in 1:tsnetwork.numcust
        r = chosen_routes[p]
        res["Passenger $p"] = routedesc[r][2]
        res1[p] = routedesc[r][2]
        walking["Passenger $p"] = walking_dict[p, r]
    end

    return res, res1, walking
end


"""
    computeArrivalTime(tsnetwork, numroutes, numstations, x, zD, Npd, R, routedesc, types, walking_speed)
    Returns Dict[p] => arrival time for passenger p
"""
function computeArrivalTime(tsnetwork, numroutes, numstations, x, zD, Npd, R, routedesc, types, walking_speed)
    res = Dict()
    chosen_routes = giveRoute(tsnetwork, x, R)
    destination_nodes = selectedNode(tsnetwork, zD, Npd) #Dict[p] => n where n is the destination node selected for passenger p if applicable (type 2, 4, 5)
    for p in 1:tsnetwork.numcust
        type = types[p]
        if type == 2 || type == 4 || type == 5
            res[p] = tsnetwork.nodedesc[destination_nodes[p]][2] #For type 2,4,5, arrival time is the time of the destination node in the OnD TSN
        end    
        if type == 1 || type == 3
            route = chosen_routes[p]
            res[p] = arrivalTimesType13(tsnetwork, numroutes, numstations, route, routedesc, walking_speed, tsnetwork.horizon, p)
        end
    end

    return res
end


"""
    computeCumulativeWaiting(tsnetwork, x, f, zO, Npo, types, R, routedesc, walking_speed, numstations)
    Returns Dict[p] => cumulative waiting time for passenger p
"""
function computeCumulativeWaiting(tsnetwork, x, f, zO, Npo, types, R, routedesc, walking_speed, numstations)
    origin_nodes = selectedNode(tsnetwork, zO, Npo)
    waiting_arcs, moving_arcs = partition_arcs(tsnetwork)
    waiting_arc_station = stationWaitingArcs(tsnetwork, waiting_arcs, numstations)
    chosen_routes = giveRoute(tsnetwork, x, R)
    walking_times_to_station = walking_t_to_station(tsnetwork, numstations, walking_speed)
    res = Dict()
    for p in 1:tsnetwork.numcust
        res[p] = 0
        type = types[p]
        r = chosen_routes[p]
        arc = routedesc[r][1]
        if arc in collect(keys(tsnetwork.transfer_arcs))
            res[p] += tsnetwork.transfer_arcs[arc][1] #Transfer waiting
        end
        if type == 2 || type == 3 || type == 5
            res[p] = tsnetwork.nodedesc[origin_nodes[p]][2] #Time before being picked up by van
        end
        for a in waiting_arc_station
            res[p] += f[p,a] #In the case van drops too soon
        end 
        if type == 1 || type == 4
            s, t = tsnetwork.nodedesc[tsnetwork.arcdesc[routedesc[r][1]][1]] 
            res[p] += (t - walking_times_to_station[p,s]) #Walking --> waiting --> Boarding
        end
    end

    return res
end


"""
    vehicleDistance(tsnetwork, y)
    Returns total distance covered by vehicles
"""
function vehicleDistance(tsnetwork, y)
    res = 0
    for a in 1:tsnetwork.numarcs
        res += y[a] * (tsnetwork.arcdesc[a][3] - 1) #to fix
    end
    return res
end

#--------TO VISUALIZE TRANSIT ARCS-------------#


"""
    inducedTransitArcs(ond_tsnetwork, x, R, routedesc)
    Returns:
    - t_arc: Set of transit arcs (indexed in the OnD TSN ! visualization purposes) that are used by at least one passenger
    - t_arc_dict: Dict[p] => Set of transit arcs (indexed in the OnD TSN ! visualization purposes) that are used by passenger p
"""
function inducedTransitArcs(ond_tsnetwork, x, R, routedesc)
    t_arc_dict = Dict()
    t_arc = []
    chosen_routes = giveRoute(ond_tsnetwork, x, R)
    for p in 1:ond_tsnetwork.numcust
        t_arc_dict[p]  = []
        r = chosen_routes[p]
        if r != 0
            arc = routedesc[r][1]
            push!(t_arc, arc)
            push!(t_arc_dict[p], arc)
        end
    end
    return t_arc, t_arc_dict
end




