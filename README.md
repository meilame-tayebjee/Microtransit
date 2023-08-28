# Microtransit

Project realized during a research internship at MIT's Operations Research Center between April and July 2023.

"Internship_report.pdf" contains all the necessary details and explanations (context, variables, model...). All the results presented in it are replicable with this code.

This code enables to:

- generate a certain number of customers, their origins and destinations, following a certain setting (Linear Expansion or Cross)
- generate an associated time-space network
- run a Mixed-Integer Linear Programming optimization model that enables to route vehicles and passengers through this time-space network, assigning each customer to a certain type of ride (transit, on-demand, hybrid)
- visualize results

Moreover, a Benders decomposition algorithm - with and without Pareto-optimal cuts - is provided.

This project is mainly coded in Julia - except for the generation of the settings that is coded in Python 3.
Besides the classic packages that can be directly installed with Julia's built-in package manager Pkg, all the JuMP models in the code use the Gurobi optimizer (https://www.gurobi.com/) and thus need a license to run.

Thereafter, we sift through all the files.

### Scripts folder

Contains everything to generate and visualize both on-demand and transit time-space networks.

### utils.jl

Contains numerous methods that deal with time-space networks: distances, nodes, arcs...

- determine_distance(tsnetwork, l1, l2)
- distance_to_station(tsnetwork, numstations)
- walking_t_to_station(tsnetwork, numstations, walking_speed)
- distance_alight_to_dest(tsnetwork, numstations)
- distance_to_dest(tsnetwork)
- transit_arclist_transitTSN(transit_tsnetwork)
- transit_times(transit_tsnetwork, num_stations)
- transit_arclist_OnD(ond_tsnetwork, transit_tsnetwork)
- partition_arcs(tsnetwork)
- stationWaitingArcs(tsnetwork, waiting_arcs, numstations)
- arc_origin_times(tsnetwork)
- arc_arrival_times(tsnetwork)
- flowBalanceNodes(tsnetwork)
- define_zones(ond_tsnetwork)
- zone(a, tsnetwork)

### preprocessing.jl

Contains the preprocessing algorithm and all the methods needed for this algorithm and all the method related to _routes_, namely:

- OnD_origins(tsnetwork, max_waiting, departure_times)
- OnD_destinations(tsnetwork, max_arrival_times, departure_times)
- buildRoutes(transit_arc_list)
- determineRouteWalkingDistance(tsnetwork, numroutes, routedesc, times_to_station, dist_alight_dest)
- determineBoardingTime(tsnetwork, arc_id)
- maxArrivalTime(tsnetwork, walking_speed)
- boardingNodeRoute(tsnetwork, numroutes, routedesc)
- reverseBoardingNodeRoute(boarding_node_route)
- alightNodeRoute(tsnetwork, numroutes, routedesc)
- reverseAlightNodeRoute(alight_node_route)
- boardingStationRoute(tsnetwork, numroutes, routedesc)
- alightingStationRoute(tsnetwork, numroutes, routedesc)
- arrivalTimesType13(tsnetwork, numroutes, numstations, route, routedesc, walking_speed, horizon, p)
- route_boarding_alighting_time(tsnetwork, boarding_node_route, alight_node_route, numroutes)

### model.jl

Contains the MILP model as well as all the methods that concern the analysis of **optimized** results:

- f_to_arctype(f, types, tsnetwork, numstations)
- y_to_arclist(v, tsnetwork)
- path(p, f, types, ond_tsnetwork, numstations)
- selectedNode(tsnetwork, z, N)
- giveRoute(tsnetwork, x, R)
- give_type_and_walk(tsnetwork, x, routedesc, R, walking_dict)
- computeArrivalTime(tsnetwork, numroutes, numstations, x, zD, Npd, R, routedesc, types, walking_speed)
- computeCumulativeWaiting(tsnetwork, x, f, zO, Npo, types, R, routedesc, walking_speed, numstations)
- vehicleDistance(tsnetwork, y)
- inducedTransitArcs(ond_tsnetwork, x, R, routedesc)

### run_optimization.jl

The main file that enables to run directly the Gurobi optimization model. Parameters can be changed in the beginning of the script.
The file calls almost all the other files of the project (except Benders.jl).
Uncomment the last line to make the run ! 

### Benders.jl 

Contains everything related to the Benders decomposition: primal and dual of the subproblem, methods to retrieve duals, the Benders algorithm itself (with and without Pareto-optimal cuts).
File is ready to be runned, but please make the last line of run_optimization.jl is commented to avoid running the latter for nothing.

### Other

generate_cross.py and generate_linear.py are the Python scripts to generate settings.
The folder visualizations contains both visualizations of the resulting TSN and the setting ; variables saves the results of the optimization so as to be able to make analysis without re-running ; and Benders saves the results of the Master Problem, the subproblem... at each iteration. All the savings are done automatically and already coded.


Please feel free to reach out to meilame.tayebjee@polytechnique.edu for questions, help or suggestions !



