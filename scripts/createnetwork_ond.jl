
#Reads the location file and returns the location coordinates
function ond_readlocations(locationfilename, numlocs)
	data = CSV.read(locationfilename, DataFrame)
	loccoords = hcat(data[:,2], data[:,3])[1:numlocs,:]
	return loccoords

end

#-----------------------------------------------------------------------------------#

#Create the time-space nodes, returning the number, ids, and descriptions
function ond_createtimespacenodes(numlocs, numcust, horizon, tstep)

	nodeid, nodedesc = Dict(), Dict()

	index = 1
	for t in 0:tstep:horizon
		for l in 1:numlocs
			nodeid[l,t] = index
			nodedesc[index] = (l,t)
			index += 1
		end	
	end
	
	numnodes = length(nodeid)
	times = [t for t in 0:tstep:horizon]

	return numnodes, nodeid, nodedesc, times

end

#-----------------------------------------------------------------------------------#

#Read the list of arcs from the arc file
function ond_getphysicalarcs(arcfilename, tstep, numlocs)

	data = CSV.read(arcfilename, DataFrame)
	
	physicalarcs = []

	for a in 1:size(data)[1]
		l1, l2 = data[a,1], data[a,2]
		if (l1 <= numlocs) & (l2 <= numlocs)
			arcdistance = data[a,4]
			arclength_raw = data[a,3]
			arclength_discretized = tstep * ceil(arclength_raw / tstep)
			push!(physicalarcs, (l1, l2, arcdistance, arclength_raw, arclength_discretized))
		end
	end

	return physicalarcs

end

#-----------------------------------------------------------------------------------#

#Create the time-space network arcs, returning the number, ids, descriptions, and cost of each arc
function ond_createtimespacearcs(physicalarcs, transit_tsnetwork, numcust, numnodes, nodeid, numlocs, horizon, tstep)

	arcid, arcdesc, A_plus, A_minus, arccost = Dict(), Dict(), Dict(), Dict(), []
	
	for node in 1:numnodes
		A_plus[node] = Set()
		A_minus[node] = Set()
	end

	stationaryarcs = []
	for l in 1:numlocs
		push!(stationaryarcs, (l, l, 0, tstep, tstep))
	end

	index = 1
	for arc in union(physicalarcs, stationaryarcs), t in 0:tstep:horizon-arc[5]
		startnode = nodeid[arc[1],t]
		endnode = nodeid[arc[2],t+arc[5]]
	
		arcid[(startnode,endnode)] = index
		arcdesc[index] = (startnode,endnode, arc[5])
		push!(A_plus[startnode], index)
		push!(A_minus[endnode], index)
		push!(arccost, arc[3])

		index += 1
	end

	numarcs = length(arcid)
	

	#Adding transit transfer arcs
	transfer_arcs = Dict()
	for arc in collect(keys(transit_tsnetwork.transfer_arcs))
		n1_transit, n2_transit, dist = transit_tsnetwork.arcdesc[arc]
        s1, time1 = transit_tsnetwork.nodedesc[n1_transit] #s1 must be in 1:num_stations
        s2, time2 = transit_tsnetwork.nodedesc[n2_transit]
        n1_ond = nodeid[(2*numcust + s1, time1)]
        n2_ond = nodeid[(2*numcust + s2, time2)]

		if (n1_ond, n2_ond) in keys(arcid)
			transfer_arcs[arcid[(n1_ond,n2_ond)]] = transit_tsnetwork.transfer_arcs[arc][1]
			continue
		end
		arcid[(n1_ond,n2_ond)] = index
		arcdesc[index] = (n1_ond, n2_ond, dist)
		push!(A_plus[n1_ond], index)
		push!(A_minus[n2_ond], index)
		push!(arccost, dist)
		transfer_arcs[index] = transit_tsnetwork.transfer_arcs[arc][1]
		index += 1
	end

	return numarcs, arcid, arcdesc, A_plus, A_minus, arccost, transfer_arcs

end

#-----------------------------------------------------------------------------------#

#Build the full time-space network
function ond_createfullnetwork(locationfilename, arcfilename, transit_tsn, numlocs, numcust, numzones, horizon, tstep)

	#Build network
	loccoords = ond_readlocations(locationfilename, numlocs)
	numnodes, nodeid, nodedesc, times = ond_createtimespacenodes(numlocs, numcust, horizon, tstep)
	physicalarcs = ond_getphysicalarcs(arcfilename, tstep, numlocs)
	numarcs, arcid, arcdesc, A_plus, A_minus, arccost, transfer_arcs = ond_createtimespacearcs(physicalarcs, transit_tsn, numcust, numnodes, nodeid, numlocs, horizon, tstep)

	#Create a NamedTuple with all the useful network data/parameters
	tsnetwork = (loccoords=loccoords, numnodes=numnodes, nodeid=nodeid, nodedesc=nodedesc, times=times, numarcs=numarcs, arcid=arcid, arcdesc=arcdesc, A_plus=A_plus, A_minus=A_minus, arccost=arccost, numlocs = numlocs, numcust = numcust, horizon = horizon, numzones = numzones, transfer_arcs = transfer_arcs, tstep = tstep)
	println("Initialized On-Demand time space network with...")
	println("Num nodes = ", tsnetwork.numnodes)
	println("Num arcs = ", tsnetwork.numarcs)
	println("Added $(length(arcid) - numarcs) transfer arcs to OnD TSN")
	return tsnetwork

end
