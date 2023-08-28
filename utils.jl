
function determine_distance(tsnetwork, l1, l2)
    #l1, l2 locations
    result = 0
    for arc in tsnetwork.A_plus[l1]
        n1, n2, dist = tsnetwork.arcdesc[arc]
        l_fin, time2 = tsnetwork.nodedesc[n2]
        if l_fin == l2
            result = dist
            break
        end
    end
    return floor(result)
end

function distance_to_station(tsnetwork, numstations)
    #Returns Dict[p, s] --> distance from p to station s loc
    result = Dict()
    for p in 1:tsnetwork.numcust
        for station in 2*tsnetwork.numcust+1:2*tsnetwork.numcust+numstations
            t_to_station = determine_distance(tsnetwork, p, station)
            result[p, station] = t_to_station
        end
    end
    return result
end

function walking_t_to_station(tsnetwork, numstations, walking_speed)
    #Returns Dict[p, s] --> time to walk from p to station s
    result = Dict()
    temp = distance_to_station(tsnetwork, numstations)
    for p in 1:tsnetwork.numcust
        for s in 2*tsnetwork.numcust+1:2*tsnetwork.numcust+numstations
            t_to_station = temp[p, s]
            result[p, s] = t_to_station / walking_speed
        end
    end
    return result
end

function distance_alight_to_dest(tsnetwork, numstations)
#Returns Dict[p, s] --> distance from s to dest(p)
    result = Dict()
    for dest in tsnetwork.numcust+1:2*tsnetwork.numcust
        for station in 2*tsnetwork.numcust+1:2*tsnetwork.numcust+numstations
            dist = determine_distance(tsnetwork, station, dest)
            result[dest - tsnetwork.numcust, station] = dist
        end
    end
    return result
end

function distance_to_dest(tsnetwork)
#Returns direct riding distance from origin to dest
    result = Dict()
    for p in 1:tsnetwork.numcust
        t_to_dest = determine_distance(tsnetwork, p, p + tsnetwork.numcust)
        result[p] = t_to_dest
    end
    return result
end


#These methods are focused on transit TSN
#----------------------------------------#
function transit_arclist_transitTSN(transit_tsnetwork)
    transit_arcs, waiting_arcs = partition_arcs(transit_tsnetwork)
    return transit_arcs
end

function transit_times(transit_tsnetwork, num_stations)
    res = Dict() # key = departure station num, value = Set of times of departure
    for i in 1:num_stations
        res[i] = Set()
    end
    transit_arcs = transit_arclist_transitTSN(transit_tsnetwork)
    for arc in transit_arcs
        n1, n2, dist = transit_tsnetwork.arcdesc[arc]
        l1, time1 = transit_tsnetwork.nodedesc[n1] #l1 must be in 1:num_stations
        push!(res[l1], time1)
    end
    return res
end
#----------------------------------------#


"""
    transit_arclist_OnD(ond_tsnetwork, transit_tsnetwork)
    This method returns the list of arcs in the OnD TSN corresponding to the transit arcs in the transit TSN, as well as the boarding nodes.
    It sort of connects both time-space networks

"""
function transit_arclist_OnD(ond_tsnetwork, transit_tsnetwork)
    arclist = Set()
    boarding_nodes_OnD = Set()
    transit_arcs_TSN = transit_arclist_transitTSN(transit_tsnetwork)
    for arc in transit_arcs_TSN
        n1_transit, n2_transit, dist = transit_tsnetwork.arcdesc[arc]
        s1, time1 = transit_tsnetwork.nodedesc[n1_transit] #s1 must be in 1:num_stations
        s2, time2 = transit_tsnetwork.nodedesc[n2_transit]
        n1_ond = ond_tsnetwork.nodeid[(2*ond_tsnetwork.numcust + s1, time1)] #2*ond_tsnetwork.numcust + s1 corresponds to the location of the same station but in the OnD TSN!
        n2_ond = ond_tsnetwork.nodeid[(2*ond_tsnetwork.numcust + s2, time2)]
        arc_ond = ond_tsnetwork.arcid[n1_ond, n2_ond]
        push!(arclist, arc_ond)
        push!(boarding_nodes_OnD, n1_ond)
    end
    return arclist, boarding_nodes_OnD
end

#Splits the arcs of the OnD TSN into 2 sets: moving arcs and waiting arcs
function partition_arcs(tsnetwork)
    moving = Set()
    waiting = Set()
    for arc in 1:tsnetwork.numarcs
        arc_tuple = tsnetwork.arcdesc[arc]
        n1, n2, dist = arc_tuple
        l1, time1 = tsnetwork.nodedesc[n1]
        l2, time2 = tsnetwork.nodedesc[n2]

        if l1 != l2
            push!(moving, arc)
        else 
            push!(waiting, arc)
        end
    end
    return moving, waiting
end

#Returns the arcs of the OnD TSN that are waiting arcs AT STATIONS
function stationWaitingArcs(tsnetwork, waiting_arcs, numstations)
    res = Set()
    stations = [2*tsnetwork.numcust+1:2*tsnetwork.numcust+numstations]
    for arc in waiting_arcs
        n1, n2, dist = tsnetwork.arcdesc[arc]
        l1, time1 = tsnetwork.nodedesc[n1]
        if l1 in stations
            push!(res, arc)
        end
    end
    return res
end

#For an arc ((l1, t1), (l2, t2)) compiles t1
#Returns Dict[arc] --> t1
function arc_origin_times(tsnetwork)
    times = Dict()
    for arc in 1:tsnetwork.numarcs
        n1, n2, dist = tsnetwork.arcdesc[arc]
        l1, time1 = tsnetwork.nodedesc[n1]
        times[arc] = time1
    end
    return times
end

#For an arc ((l1, t1), (l2, t2)) compiles t2
#Returns Dict[arc] --> t2
function arc_arrival_times(tsnetwork)
    times = Dict()
    for arc in 1:length(tsnetwork.arcid)              
        n1, n2, dist = tsnetwork.arcdesc[arc]
        l2, time2 = tsnetwork.nodedesc[n2]
        times[arc] = time2
    end
    return times
end

#It is the set \mathcal{N}_0 in the paper
function flowBalanceNodes(tsnetwork)
    result = Set()
    for i in 1:tsnetwork.numlocs
        for t in 1:tsnetwork.horizon-1
            push!(result, tsnetwork.nodeid[i, t])
        end
    end

    return result
end

#Here: 1 zone = 1 station
function define_zones(ond_tsnetwork)
    I = Dict()
    for i in 1:ond_tsnetwork.numzones
        I[i] = []
        push!(I[i], 2*ond_tsnetwork.numcust + i)

    end
    return I
end


function typesRepartitionToString(types_repartition::Dict{Any, Any})
    sortedKeys = sort(collect(keys(types_repartition)), rev = false)
    repartition = join([string(round(types_repartition[key], digits=2)) for key in sortedKeys], "/")
    return repartition
end

"""
    indicator(a, ensemble)
    Returns 1 if a ∈ ensemble (a can be anything)

TBW
"""
function indicator(a, ensemble)
    if a in ensemble
        return 1
    else
        return 0
    end
end


"""
    zone(a, tsnetwork)
    Given an arc ((l1, t1), (l2, t2)), this method returns the zone of location l1 in the OnD TSN
"""
function zone(a, tsnetwork)
    n1, n2, dist = tsnetwork.arcdesc[a]
    l1, time1 = tsnetwork.nodedesc[n1]
    return l1 - 2*tsnetwork.numcust
end

