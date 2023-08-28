using Luxor, Colors

function matchingTransitTSNtoOnDTSN(transit_tsnetwork, ond_tsnetwork)
    res = Dict()
    for i in 1:transit_tsnetwork.numarcs
        n1, n2, dist = transit_tsnetwork.arcdesc[i]
        s1, t1 = transit_tsnetwork.nodedesc[n1]
        s2, t2 = transit_tsnetwork.nodedesc[n2]
        n1_ond = ond_tsnetwork.nodeid[(2*ond_tsnetwork.numcust + s1, t1)]
        n2_ond = ond_tsnetwork.nodeid[(2*ond_tsnetwork.numcust + s2, t2)]
        res[i] = ond_tsnetwork.arcid[n1_ond, n2_ond]
    end
    return res
end

function matchingOnDTSNtoTransitTSN(transit_tsnetwork, ond_tsnetwork, transit_arclist_ond)
    res = Dict()
    for i in transit_arclist_ond
        n1, n2, dist = ond_tsnetwork.arcdesc[i]
        l1, t1 = ond_tsnetwork.nodedesc[n1]
        l2, t2 = ond_tsnetwork.nodedesc[n2]
		if l1 > 2*ond_tsnetwork.numcust && l2 > 2*ond_tsnetwork.numcust
			n1_trans = transit_tsnetwork.nodeid[(l1 - 2*ond_tsnetwork.numcust, t1)]
			n2_trans = transit_tsnetwork.nodeid[(l2 - 2*ond_tsnetwork.numcust, t2)]
			res[i] = transit_tsnetwork.arcid[n1_trans, n2_trans]
		end
    end
    return res
end


function typeColorMatching()
	res = Dict()
	res[1] = (0, 0, 0) #Should not appear. If it does, MOVING_COST is most likely too low
	res[2] = (155, 0, 0)
	res[3] = (100, 125, 255)
	res[4] = (50, 255, 175)
	res[5] = (100, 100, 50)

	return res
end

function ond_timespaceviz(drawingname, tsn, transit_tsn, arc_to_type, arclist2, transitarcs; x_size=1200, y_size=700)
	#transit arcs are the transit arcs in the OnD TSN. See method transit_arclist_OnD in utils.jl
	matching_OnD_transit = matchingOnDTSNtoTransitTSN(transit_tsn, tsn, transitarcs)
	matching_transit_OnD = matchingTransitTSNtoOnDTSN(transit_tsn, tsn)

	numcust, numlocs = tsn.numcust, tsn.numlocs
	#Find coordinates for each time-space node
	nodelist = []
	x_size_trimmed, y_size_trimmed = x_size*0.9, y_size*0.9
	k1 = x_size_trimmed/(tsn.horizon/tsn.tstep + 2) 
	k2 = y_size_trimmed/(tsn.numlocs + 2)
	for i in 1:tsn.numnodes
		ycoord = tsn.nodedesc[i][1]
		xcoord = (tsn.nodedesc[i][2]/tsn.tstep)+1

		#Scaling to image size
		tup = (-x_size_trimmed/2 + xcoord*k1, -y_size_trimmed/2 + ycoord*k2)   
		
		push!(nodelist,tup)
	end

	#Create actual points as a Luxor object
	nodePoints = Point.(nodelist)

	#---------------------------------------------------------------------------------------#

	#Arcs for visualization
	#Duplicate for multiple input arc lists with different colors/thickness/dash if you're trying to show m
	matching = typeColorMatching()
	arclist = collect(keys(arc_to_type))
	arcinfo = []
	for a in arclist
		startPoint = nodePoints[tsn.arcdesc[a][1]]
		endPoint = nodePoints[tsn.arcdesc[a][2]]
		
		#Set arc attributes
		arcColor = matching[arc_to_type[a]] #RGB tuple 
		arcDash = "solid" #"solid", "dashed"			
		arcThickness = 4 
		
		#Add to arcinfo list to be used in the drawing 
		push!(arcinfo, (startPoint, endPoint, arcColor, arcDash, arcThickness))
	end

	for a in arclist2
		startPoint = nodePoints[tsn.arcdesc[a][1]]
		endPoint = nodePoints[tsn.arcdesc[a][2]]
		
		#Set arc attributes
		arcColor = (255,0,255) #RGB tuple 
		arcDash = "dashed" #"solid", "dashed"			
		arcThickness = 3
		
		#Add to arcinfo list to be used in the drawing 
		push!(arcinfo, (startPoint, endPoint, arcColor, arcDash, arcThickness))
	end

	for a in transitarcs
		transit_a = matching_OnD_transit[a]
		#Set arc attributes
		if transit_a in collect(keys(transit_tsn.transfer_arcs))
				for transit_arc_temp in transit_tsn.transfer_arcs[transit_a][2]
					arc_temp = matching_transit_OnD[transit_arc_temp]
					startPoint = nodePoints[tsn.arcdesc[arc_temp][1]]
					endPoint = nodePoints[tsn.arcdesc[arc_temp][2]]
					arcColor = (255,255,50) #RGB tuple 
					arcDash = "solid" #"solid", "dashed"			
					arcThickness = 4 
					push!(arcinfo, (startPoint, endPoint, arcColor, arcDash, arcThickness))
				end
			else
				#Set arc attributes
				startPoint = nodePoints[tsn.arcdesc[a][1]]
				endPoint = nodePoints[tsn.arcdesc[a][2]]
				arcColor = (0,255,0) #RGB tuple 
				arcDash = "solid" #"solid", "dashed"			
				arcThickness = 4 
				push!(arcinfo, (startPoint, endPoint, arcColor, arcDash, arcThickness))
			end

		#Add to arcinfo list to be used in the drawing 
		push!(arcinfo, (startPoint, endPoint, arcColor, arcDash, arcThickness))
	end

	#-------------------------------------------------------------------------#

	#Initiailize drawing
	Drawing(x_size, y_size, drawingname)
	origin()
	background("white")

	#Draw arcs
	for i in arcinfo
		
		#Set arc attributes from the arcinfo
		r_val, g_val, b_val = i[3][1]/255, i[3][2]/255, i[3][3]/255
		setcolor(convert(Colors.HSV, Colors.RGB(r_val, g_val, b_val)))  #You can also use setcolor("colorname")
		setdash(i[4])
		setline(i[5])

		#Draw the line from the start node to end node
		line(i[1], i[2] , :stroke)
		
		#Figure out the angle of the arrow head
		theta = atan((i[2][2] - i[1][2])/(i[2][1] - i[1][1]))
		dist = distance(i[1], i[2])
		arrowhead = (1-8/dist)*i[2] + (8/dist)*i[1] #8 pixels from the end node
		
		#Draw the arrow head
		local p = ngon(arrowhead, 5, 3, theta, vertices=true)
		poly(p, :fill,  close=true)
	end

	#Draw node points
	setcolor("black")
	circle.(nodePoints, 4, :fill)

	#Set font size for labels
	fontsize(14)

	#Add location labels
	for l in 1:tsn.numcust
		coord = nodePoints[tsn.nodeid[(l,0.0)]]
		label("Origin $l       ", :W , coord)
	end
	for l in numcust+1:2*numcust
		coord = nodePoints[tsn.nodeid[(l,0.0)]]
		o = l - numcust
		label("Destination $o   ", :W , coord)
	end
	for l in 2*numcust+1:numlocs
		o = l - 2*numcust
		coord = nodePoints[tsn.nodeid[(l,0.0)]]
		label("Station $o ", :W , coord)
	end

	#Add time labels
	for t in 0:tsn.tstep*2:tsn.horizon
		coord = nodePoints[tsn.nodeid[(1,t)]] + Point(0,-30)
		label("t = $t", :N , coord)
	end

	finish()
	preview()

end
