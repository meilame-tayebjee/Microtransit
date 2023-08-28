using  Random, CSV, DataFrames, Statistics, Dates, PyCall, JLD
randomseedval = 1			#Set random seed if there are any random components of an algorithm
Random.seed!(randomseedval)

#-------------------------------------LOAD SCRIPTS-------------------------------------#

include("scripts/createnetwork_ond.jl")
include("scripts/networkvisualization_ond.jl")
include("scripts/createnetwork_transit.jl")
include("scripts/networkvisualization_transit.jl")
include("model.jl")

#-------------------------------------CHOOSE SETTING-------------------------------------#

setting = "CROSS" #"CROSS" or "LINEAR_EXPANSION"

#----------------------------------PARAMETERS----------------------------------#  	

SAMPLE_SIZE = 1
#Read experiment parameters 
horizon = 100									#Length of time horizon in hours 
tstep = 1       								#Time discretization

numcust = 0
if setting == "CROSS"
    numcust = 4*SAMPLE_SIZE                      #Number of customers
    numstations = 5                
elseif setting == "LINEAR_EXPANSION"
    numcust = SAMPLE_SIZE
    numstations = 3                    
end

numlocs = 2*numcust + numstations				#Number of physical locations
numzones = numstations

ond_locationfilename = "data/locations_$SAMPLE_SIZE.csv"
ond_arcfilename = "data/arcs_$SAMPLE_SIZE.csv"
transit_locationfilename = "data/transit_locations.csv"
transit_arcfilename = "data/transit_arcs.csv"

# PREPROCESSING PARAMETERS #
MAX_WAITING_TRANSFER = 2       #Max waiting time allowed for a transfer
WALKING_SPEED = 0.5             #To be selected within ]0, 1[ as we walk slower than the vans...
MAX_WAITING = 20                #Max cumulated waiting allowed for a customer
MAX_WALKING = 50                 #Max cumulated walking allowed for a customer

# TRANSPORTATION PARAMETERS #
freq =10                         #Frequency of transit lines
WALKING_COST = 80               #"Pain" of walking per km
MOVING_COST =  10               #Cost per km for on-demand vehicles
K = 4                           #Vans capacity
W = [3, 3, 3, 3, 0]            #Number of vans PER ZONE ; array size = numzones
if length(W) != numzones
    println("Error: W array size is not equal to numzones")
end
cap = sum(W)                    #Total number of vans
epsilon = 0.5                   #Steady state parameter

TIME_COSTS = Dict()
for t in 0:horizon
    TIME_COSTS[t] = sqrt(t)           #"Time pressure". Choose between sqrt(t), log(t) and t. /!\ MOVING_COST might have to be updated to avoid absurd results
end


DEPARTURE_INTERVAL = 5      #Size of the departure window
nb_iter = 1                 #Monte-Carlo iterations
deltas = [20]               #Deltas we want to test
radiuses = [5]              #Radiuses we want to test

#-------------------------------------RUN SCRIPTS-------------------------------------#
function run_opti(nb_iter, deltas, radiuses)


    departure_times = rand(1:DEPARTURE_INTERVAL, numcust)
    println("Departure times = ", departure_times)
    df = DataFrame(Delta = Int[], Radius = Int[], AverageWaiting = Float64[], AverageWalking = Float64[], AverageArrival = Float64[], TypesRepartition = String[], PreprAvgTime = Float64[], OptAvgTime = Float64[])


    for delta in deltas
        for radius in radiuses 
            RADIUS = radius
            DELTA = delta #Length of transit line
            DELTA2 = DELTA


            waitings = []
            walkings = []
            arrivals = []
            preprocessing_times = []
            optimizing_times = []
            types_repartition = Dict()
            for i in 1:5
                types_repartition[i] = 0
            end
            try
                for iter in 1:nb_iter
                    #-------------------------------------GENERATE SETTING-------------------------------------#
                    #Replace "python" by "python3" if needed
                    if setting == "CROSS"
                        run(`python generate_cross.py $RADIUS $SAMPLE_SIZE $DELTA $DELTA2 $iter`)
                    elseif setting == "LINEAR_EXPANSION"
                        run(`python generate_linear.py $RADIUS $SAMPLE_SIZE $DELTA $DELTA2 $iter`)
                    end

                    #-----------------------------------GENERATE NETWORKS-----------------------------------# 
                    #Create node and arc networks
                    transit_tsnetwork = transit_createfullnetwork(transit_locationfilename, transit_arcfilename, numstations, horizon, tstep, freq, MAX_WAITING_TRANSFER, setting)
                    ond_tsnetwork = ond_createfullnetwork(ond_locationfilename, ond_arcfilename, transit_tsnetwork, numlocs, numcust, numzones, horizon, tstep)

                    #----------------------------------RUN PREPROCESSING----------------------------------#
                    #Routes preprocessing
                    #To be noted that R2 contains only route n° 0 which encodes the full on-demand route (no transit arc)
                    R, R1, R2, R3, R4, R5, Npo, Npb, Npa, Npd, routeid, routedesc, numroutes, walking_dict = preprocessing(ond_tsnetwork, transit_tsnetwork, numstations, MAX_WALKING, MAX_WAITING, WALKING_SPEED, departure_times)

                    #--------------------------------OPTIMIZATION---------------------------------#
            
                    f, y , x , zO, zD, preprocessing_t, optimizing_t = model(ond_tsnetwork, numcust, numstations, numzones, horizon,  K, W, R, R1, R2, R3, R4, R5, Npo, Npb, Npa, Npd, WALKING_SPEED,routedesc, numroutes, TIME_COSTS, MOVING_COST, WALKING_COST)

                    types_p, types, walking = give_type_and_walk(ond_tsnetwork, x, routedesc, R, walking_dict)
                    save("variables/variables_$iter.jld", "f", f, "y", y, "x", x, "zO", zO, "zD", zD)

                    #--------------------------------ANALYSIS---------------------------------# 

                    arrival_times = computeArrivalTime(ond_tsnetwork, numroutes, numstations, x, zD, Npd, R, routedesc, types, WALKING_SPEED)
                    waiting = computeCumulativeWaiting(ond_tsnetwork, x, f, zO, Npo,types, R, routedesc, WALKING_SPEED, numstations)
                    distance_veh = vehicleDistance(ond_tsnetwork, y)
                    println("")
                    println("------ SOME STATISTICS ------")
                    println("Types = ", types_p)
                    println("")
                    println("Walking distance = ", walking)
                    println("Total walking distance = ", sum(collect(values(walking))))
                    println("")
                    println("Arrival times = ", arrival_times)
                    println("")
                    println("Waiting times = ", waiting)
                    println("Total waiting time = ", sum(collect(values(waiting))))
                    println("")
                    println("Distance total covered by vehicles = ", distance_veh )

                    #---------------------------------MONTE-CARLO COMPUTATION---------------------------------#
                    push!(waitings, sum(collect(values(waiting))))
                    push!(walkings, sum(collect(values(walking))))
                    push!(arrivals, maximum(collect(values(arrival_times))))
                    push!(preprocessing_times, preprocessing_t)
                    push!(optimizing_times, optimizing_t)
                    for p in 1:numcust
                        types_repartition[types[p]] += 1/nb_iter
                    end


                    #--------------------------------VISUALIZATION---------------------------------# 

                    arc_to_type = f_to_arctype(f, types, ond_tsnetwork, numstations)
                    arclist2 = y_to_arclist(y, ond_tsnetwork)
                    t_arc, t_arc_dict = inducedTransitArcs(ond_tsnetwork, x, R, routedesc)

                    filename = "visualizations/ondemand_tsn" * "$iter.png"
                    println("ITERATION $iter DONE")
                    ond_timespaceviz(filename, ond_tsnetwork, transit_tsnetwork, arc_to_type , arclist2, t_arc, x_size=2000, y_size=1000)


                end

                push!(df, (delta, radius, mean(waitings), mean(walkings), mean(arrivals), typesRepartitionToString(types_repartition), mean(preprocessing_times), mean(optimizing_times)))
                println(df)
                CSV.write("size_$SAMPLE_SIZE.csv", df)

                catch e
                    println("ERROR : ", e)
                    println("Delta ", delta, " Radius ", radius)
                    push!(df, (delta, radius, 0, 0, 0, "ERROR", 0, 0))
                    println(df)
                    CSV.write("size_$SAMPLE_SIZE.csv", df)
                    continue
                end 
        end
    end
    return departure_times
end

#UNCOMMENT NEXT LINE AND RUN TO LAUNCH MILP OPTIMIZATION 
#departure_times = run_opti(nb_iter, deltas, radiuses)
