#Reads the location file and returns the location coordinates
function readlocations(locationfilename, numstations)
	data = CSV.read(locationfilename, DataFrame)
	loccoords = hcat(data[:,2], data[:,3])[1:numstations,:]
	return loccoords

end

#-----------------------------------------------------------------------------------#

#Create the time-space nodes, returning the number, ids, and descriptions
function createtimespacenodes(numstations, horizon, tstep)

	nodeid, nodedesc = Dict(), Dict()

	index = 1
	for t in 0:tstep:horizon
		for l in 1:numstations
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
function getphysicalarcs(arcfilename, tstep, numstations, setting)

	data = CSV.read(arcfilename, DataFrame)
	
	physicalarcs = []

	for a in 1:size(data)[1]
		l1, l2 = data[a,1], data[a,2]
		if (l1 <= numstations) & (l2 <= numstations)
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
function transit_createtimespacearcs(horizon, tstep, physicalarcs, numstations, numnodes, nodeid, nodedesc, freq, max_waiting_transfer)

	arcid, arcdesc, A_plus, A_minus, arccost = Dict(), Dict(), Dict(), Dict(), []
	for node in 1:numnodes
		A_plus[node] = []
		A_minus[node] = []
	end

	stationaryarcs = []
	for l in 1:numstations
		push!(stationaryarcs, (l, l, 0, tstep, tstep))
	end

	index = 1

	for arc in stationaryarcs, t in 0:tstep:horizon-arc[5]
		startnode = nodeid[arc[1],t]
		endnode = nodeid[arc[2],t+arc[5]]
	
		arcid[(startnode,endnode)] = index
		arcdesc[index] = (startnode,endnode, arc[5])
		push!(A_plus[startnode], index)
		push!(A_minus[endnode], index)

		push!(arccost, arc[3])

		index += 1	
	end

	if setting == "CROSS"
		counter = 0
		for transitarc in union(physicalarcs[1:1:2], physicalarcs[7:1:8])
			DELTA = transitarc[5]
			for t in freq:freq:horizon-(1+counter)*transitarc[5]
				startnode = nodeid[transitarc[1], Int(t + counter * DELTA)]
				endnode = nodeid[transitarc[2],t+transitarc[5] + counter * DELTA]
			
				arcid[(startnode,endnode)] = index
				arcdesc[index] = (startnode,endnode, transitarc[5])
				push!(A_plus[startnode], index)
				push!(A_minus[endnode], index)

				push!(arccost, transitarc[3])

				index += 1	
			end
		end

		counter = 0
		for transitarc in physicalarcs[3:1:4]
			DELTA = transitarc[5]
			for t in freq:freq:horizon-(1+counter)*transitarc[5]
				startnode = nodeid[transitarc[1], Int(t + counter * DELTA)]
				endnode = nodeid[transitarc[2],t+transitarc[5] + counter * DELTA]
				arcid[(startnode,endnode)] = index
				arcdesc[index] = (startnode,endnode, transitarc[5])
				push!(A_plus[startnode], index)
				push!(A_minus[endnode], index)

				push!(arccost, transitarc[3])

				index += 1	
			end
			counter += 1
		end
		counter = 0
		for transitarc in physicalarcs[9:1:10]
			DELTA = transitarc[5]
			for t in freq:freq:horizon-(1+counter)*transitarc[5]
				startnode = nodeid[transitarc[1], Int(t + counter * DELTA)]
				endnode = nodeid[transitarc[2],t+transitarc[5] + counter * DELTA]
			
				arcid[(startnode,endnode)] = index
				arcdesc[index] = (startnode,endnode, transitarc[5])
				push!(A_plus[startnode], index)
				push!(A_minus[endnode], index)

				push!(arccost, transitarc[3])

				index += 1	
			end
			counter += 1
		end

		counter = 0
		for transitarc in physicalarcs[5:1:6]
			DELTA = transitarc[5]
			for t in freq:freq:horizon-(1+counter)*transitarc[5]
				startnode = nodeid[transitarc[1], Int(t + counter * DELTA)]
				endnode = nodeid[transitarc[2],t+transitarc[5] + counter * DELTA]
			
				arcid[(startnode,endnode)] = index
				arcdesc[index] = (startnode,endnode, transitarc[5])
				push!(A_plus[startnode], index)
				push!(A_minus[endnode], index)

				push!(arccost, transitarc[3])

				index += 1	
			end
			counter += 1
		end
		counter = 0
		for transitarc in physicalarcs[11:1:12]
			DELTA = transitarc[5]
			for t in freq:freq:horizon-(1+counter)*transitarc[5]
				startnode = nodeid[transitarc[1], Int(t + counter * DELTA)]
				endnode = nodeid[transitarc[2],t+transitarc[5] + counter * DELTA]
			
				arcid[(startnode,endnode)] = index
				arcdesc[index] = (startnode,endnode, transitarc[5])
				push!(A_plus[startnode], index)
				push!(A_minus[endnode], index)

				push!(arccost, transitarc[3])

				index += 1	
			end
			counter += 1
		end
	elseif setting == "LINEAR_EXPANSION"
		counter = 0
		DELTA = physicalarcs[1][5]
		for transit_arcs in physicalarcs
			for t in freq:freq:horizon-(1+counter)*DELTA
				startnode = nodeid[transit_arcs[1], Int(t + counter * DELTA)]
				endnode = nodeid[transit_arcs[2], Int(t+transit_arcs[5] + counter * DELTA)]
			
				arcid[(startnode,endnode)] = index
				arcdesc[index] = (startnode,endnode, transit_arcs[5])
				push!(A_plus[startnode], index)
				push!(A_minus[endnode], index)

				push!(arccost, transit_arcs[3])

				index += 1	
			end
			counter += 1
		end
	end
	numarcs = length(arcid)

	#Adding "transfer arcs"

	connected = Dict()
	for l1 in 1:numstations
		for l2 in 1:numstations
			if l1 != l2
				connected[(l1,l2)] = false
			end
		end
	end
	for arc in 1:numarcs
		n1, n2, dist = arcdesc[arc]
		l1, t1 = nodedesc[n1]
		l2, t2 = nodedesc[n2]
		if l1 != l2
			connected[(l1,l2)] = true
		end
	end

	transfer_arcs = Dict()
	for i in 1:numarcs
		for j in 1:numarcs
			n1, n2, dist1 = arcdesc[i]
			n3, n4, dist2 = arcdesc[j]
			l1, t1 = nodedesc[n1]
			l2, t2 = nodedesc[n2]
			l3, t3 = nodedesc[n3]
			l4, t4 = nodedesc[n4]
			if 0<= t3-t2 <= max_waiting_transfer
				if !((n1, n4) in collect(keys(arcid)))
					if l1 != l2 && l2 == l3 && l3 != l4 && l1 != l4
						if connected[l1, l4] == false
							waiting_time = t3 - t2
							arcid[(n1,n4)] = index
							arcdesc[index] = (n1,n4, dist1 + dist2 + waiting_time)
							push!(A_plus[n1], index)
							push!(A_minus[n4], index)
							push!(arccost, dist1 + dist2 + waiting_time)
							intermediary_arcs = union([i,j], [arcid[(nodeid[l2, t], nodeid[l2, t + tstep])] for t in t2:tstep:t3-tstep])
							transfer_arcs[index] = (waiting_time, intermediary_arcs)
							index += 1
						end
					end
				end
			end
		end
	end
	numarcs = length(arcid)
	return numarcs, arcid, arcdesc, A_plus, A_minus, arccost, transfer_arcs
end

#-----------------------------------------------------------------------------------#

#Build the full time-space network
function transit_createfullnetwork(locationfilename, arcfilename, numstations, horizon, tstep, freq, max_waiting_transfer, setting)

	#Build network
	loccoords = readlocations(locationfilename, numstations)
	numnodes, nodeid, nodedesc, times = createtimespacenodes(numstations, horizon, tstep)
	physicalarcs = getphysicalarcs(arcfilename, tstep, numstations, setting)
	numarcs, arcid, arcdesc, A_plus, A_minus, arccost, transfer_arcs = transit_createtimespacearcs(horizon, tstep, physicalarcs, numstations, numnodes, nodeid, nodedesc, freq, max_waiting_transfer)

	#Create a NamedTuple with all the useful network data/parameters
	tsnetwork = (loccoords=loccoords, numnodes=numnodes, nodeid=nodeid, nodedesc=nodedesc, times=times, numarcs=numarcs, arcid=arcid, arcdesc=arcdesc, A_plus=A_plus, A_minus=A_minus, arccost=arccost, freq=freq, transfer_arcs = transfer_arcs)

	println("Initialized Transit time space network with...")
	println("Num nodes = ", tsnetwork.numnodes)
	println("Num arcs = ", tsnetwork.numarcs)
	println("Among which num transfer arcs = ", length(tsnetwork.transfer_arcs))

	return tsnetwork

end