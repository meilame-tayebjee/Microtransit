using JLD2, FileIO, Suppressor, Plots
include("model.jl")
include("scripts/networkvisualization_ond.jl")
include("run_optimization.jl")

randomseedval = 1			#Set random seed if there are any random components of an algorithm
Random.seed!(randomseedval)


"""
    subproblem_primal(ond_tsnetwork, numcust, numstations, K, W, R3, R4, R5, Npo, Npb, Npa, Npd, walking_speed, routedesc, numroutes, MOVING_COST, x, zO, zD, time_limit)
    The primal of the subproblem
    Returns a JuMP model
"""
function subproblem_primal(ond_tsnetwork, numcust, numstations, K, W, R3, R4, R5, Npo, Npb, Npa, Npd, walking_speed, routedesc, numroutes, MOVING_COST, x, zO, zD, time_limit)
    model = Model(Gurobi.Optimizer)
    set_optimizer_attributes(
            model, 
            "TimeLimit" => time_limit, 
            "InfUnbdInfo" => 1, 
            "DualReductions" => 0, 
            "OutputFlag" => 1,
        )

    #---------------------UTILS----------------------------#
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
    walking_times_to_station = walking_t_to_station(ond_tsnetwork, numstations, walking_speed)

    moving_arcs, waiting_arcs = partition_arcs(ond_tsnetwork)
    flow_balance_nodes = flowBalanceNodes(ond_tsnetwork) #for vehicles

    dist_alight_to_dest = distance_alight_to_dest(ond_tsnetwork, numstations) #Dict[p, s] --> distance from station s to destination of p
    station_waiting_arcs = stationWaitingArcs(ond_tsnetwork, waiting_arcs,numstations)

    f_arc_keys = union(moving_arcs, station_waiting_arcs)

    nontransfer_arcs = 1:ond_tsnetwork.numarcs

    routeWalkingDict = determineRouteWalkingDistance(ond_tsnetwork, numroutes, routedesc, walking_times_to_station, dist_alight_to_dest)
    boarding_route = boardingNodeRoute(ond_tsnetwork, numroutes, routedesc)
    alighting_route = alightNodeRoute(ond_tsnetwork, numroutes, routedesc)
    reverse_boarding_route = reverseBoardingNodeRoute(boarding_route)
    reverse_alighting_route = reverseAlightNodeRoute(alighting_route)

    #--------------------- VARIABLES -------------------------------# 

    @variable(model, 0<= y[a in nontransfer_arcs])
    @variable(model, 0 <= f[p in 1:numcust, a in f_arc_keys] <= 1)



    #--------------------- CONSTRAINTS -------------------------------#

    #Constraints 31
    @constraint(model, [p in 1:numcust, n in Npo[p]], sum(f[p,a] for a in intersect(ond_tsnetwork.A_plus[n], f_arc_keys)) - sum(f[p,a] for a in intersect(ond_tsnetwork.A_minus[n], f_arc_keys)) == zO[p,n])
    @constraint(model, [p in 1:numcust, n in Npb[p]], sum(f[p,a] for a in intersect(ond_tsnetwork.A_plus[n], f_arc_keys)) - sum(f[p,a] for a in intersect(ond_tsnetwork.A_minus[n], f_arc_keys)) == - sum(x[p,r] for r in intersect(union(R3[p], R5[p]), reverse_boarding_route[n])))
    @constraint(model, [p in 1:numcust, n in Npa[p]], sum(f[p,a] for a in intersect(ond_tsnetwork.A_plus[n], f_arc_keys)) - sum(f[p,a] for a in intersect(ond_tsnetwork.A_minus[n], f_arc_keys)) == sum(x[p,r] for r in intersect(union(R4[p], R5[p]), reverse_alighting_route[n])))
    @constraint(model, [p in 1:numcust, n in Npd[p]], sum(f[p,a] for a in intersect(ond_tsnetwork.A_plus[n], f_arc_keys)) - sum(f[p,a] for a in intersect(ond_tsnetwork.A_minus[n], f_arc_keys)) == - zD[p,n])
    @constraint(model, [p in 1:numcust, n in other_nodes[p]], sum(f[p,a] for a in intersect(ond_tsnetwork.A_plus[n], f_arc_keys)) - sum(f[p,a] for a in intersect(ond_tsnetwork.A_minus[n], f_arc_keys)) ==0)



    #Vehicle constraints
    @constraint(model, [a in f_arc_keys], sum(f[p,a] for p in 1:numcust) <= K * y[a]) #34
    @constraint(model, [n in flow_balance_nodes], sum(y[a] for a in intersect(ond_tsnetwork.A_plus[n], nontransfer_arcs)) - sum(y[a] for a in intersect(ond_tsnetwork.A_minus[n], nontransfer_arcs)) == 0) #35

    for i in 1:length(I)
        @constraint(model, [i:i], sum(sum(y[a] for a in intersect(ond_tsnetwork.A_plus[n], nontransfer_arcs)) for n in I[i]) == W[i])
    end
    #@constraint(model, [n in 1:2*ond_tsnetwork.numcust], sum(y[a] for a in intersect(ond_tsnetwork.A_plus[n], nontransfer_arcs)) == 0) #37


    @expression(model, moving_cost_expr, sum(f[p, a] * ond_tsnetwork.arccost[a] for p in 1:numcust for a in moving_arcs))
    @expression(model, nontransfer_cost_expr, sum(y[a] * ond_tsnetwork.arccost[a] for a in intersect(moving_arcs, nontransfer_arcs)))

    @objective(model, Min, MOVING_COST * (moving_cost_expr + nontransfer_cost_expr))

    return model
end

"""
    retrieve_duals(model, ond_tsnetwork, numcust, Npo, Npb, Npa, Npd, other_nodes, flow_balance_nodes, nontransfer_arcs, f_arc_keys, I)
    Once optimized, retrieve the duals of the primal subproblem
"""
function retrieve_duals(model, ond_tsnetwork, numcust, Npo, Npb, Npa, Npd, other_nodes, flow_balance_nodes, nontransfer_arcs, f_arc_keys, I)
    duals = Float64[]
    for constr in JuMP.all_constraints(model, include_variable_in_set_constraints = true)
        dual_value = JuMP.dual(constr)
        push!(duals, dual_value)
    end

    alpha_true = Array{Float64}(undef,numcust, ond_tsnetwork.numnodes)

    i = 1
    for p in 1:numcust
        for n in Npo[p]
            alpha_true[p,n] = duals[i]
            i += 1
        end
    end

    for p in 1:numcust
        for n in Npb[p]
            alpha_true[p,n] = duals[i]
            i += 1
        end
    end

    for p in 1:numcust
        for n in Npa[p]
            alpha_true[p,n] = duals[i]
            i += 1
        end
    end

    for p in 1:numcust
        for n in Npd[p]
            alpha_true[p,n] = duals[i]
            i += 1
        end
    end

    for p in 1:numcust
        for n in other_nodes[p]
            alpha_true[p,n] = duals[i]
            i += 1
        end
    end

    println("alpha done")

    lambda_true = Dict()
    for n in flow_balance_nodes
        lambda_true[n] = duals[i]
        i += 1
    end

    println("lambda done")
    gamma_true  = Array{Float64}(undef, length(I))
    for j in 1:length(I)
        gamma_true[j] = duals[i]
        i += 1
        j+=1
    end

    println("gamma done")
    beta_true = Dict()
    for a in f_arc_keys
        beta_true[a] = duals[i]
        i += 1
    end

    println("beta done")
    for a in nontransfer_arcs
        i += 1
    end

    for p in 1:numcust
        for a in f_arc_keys
            i += 1
        end
    end

    ksi_true = Dict()
    for p in 1:numcust
        for a in f_arc_keys
            ksi_true[p,a] = duals[i]
            i += 1
        end
    end
    println("ksi done")

    return alpha_true, beta_true, lambda_true, gamma_true, ksi_true
 
end


"""
    subproblem_dual(ond_tsnetwork, numcust, numzones, K, MOVING_COST, optimality_gap, moving_arcs, waiting_arcs, nontransfer_arcs, flow_balance_nodes, f_arc_keys, non_f_arc_keys, station_waiting_arcs)
    The dual of the subproblem
    Returns a JuMP model
"""
function subproblem_dual(ond_tsnetwork, numcust, numzones, K, MOVING_COST, optimality_gap, moving_arcs, waiting_arcs, nontransfer_arcs, flow_balance_nodes, f_arc_keys, non_f_arc_keys, station_waiting_arcs)

    model_dual = Model(Gurobi.Optimizer)
        @suppress set_optimizer_attributes(
                model_dual, 
                "MIPGap" => optimality_gap,
                "InfUnbdInfo" => 1, 
                "DualReductions" => 0, 
                "OutputFlag" => 0,
            )
    #---------------------UTILS----------------------------#
    t1s = arc_origin_times(ond_tsnetwork)
    gamma_arcs = Set()
    for arc in 1:ond_tsnetwork.numarcs
        if zone(arc, ond_tsnetwork) in 1:numzones 
            if t1s[arc] == 0
                push!(gamma_arcs, arc)
            end
        end
    end
    non_gamma_arcs = []
    for arc in 1:ond_tsnetwork.numarcs
        if !(arc in gamma_arcs)
            push!(non_gamma_arcs, arc)
        end
    end
    nontransfer_moving_arcs = intersect(moving_arcs, nontransfer_arcs)
    nontransfer_waiting_arcs = intersect(waiting_arcs, nontransfer_arcs)

    non_flow_balance_nodes = []
    for n in 1:ond_tsnetwork.numnodes
        if !(n in flow_balance_nodes)
            push!(non_flow_balance_nodes, n)
        end
    end

    #--------------------- VARIABLES -------------------------------# 

    @variable(model_dual, alpha[p in 1:numcust, n in 1:ond_tsnetwork.numnodes])
    @variable(model_dual, 0<= beta[a in f_arc_keys])
    @variable(model_dual, lambda[n in 1:ond_tsnetwork.numnodes])
    @variable(model_dual, gamma[k in 1:numzones])
    @variable(model_dual, 0<= ksi[p in 1:numcust, a in f_arc_keys])

    #--------------------- CONSTRAINTS -------------------------------#
    @constraint(model_dual, [n in non_flow_balance_nodes], lambda[n] == 0)

    @constraint(model_dual, [p in 1:numcust, a in moving_arcs], alpha[p, ond_tsnetwork.arcdesc[a][1]] - alpha[p, ond_tsnetwork.arcdesc[a][2]] - beta[a] - ksi[p, a] <= MOVING_COST * ond_tsnetwork.arccost[a])
    @constraint(model_dual, [p in 1:numcust, a in station_waiting_arcs], alpha[p, ond_tsnetwork.arcdesc[a][1]] - alpha[p, ond_tsnetwork.arcdesc[a][2]] - beta[a] - ksi[p, a] <= 0)

    @constraint(model_dual, [a in intersect(nontransfer_moving_arcs, f_arc_keys, gamma_arcs)], K*beta[a] + lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes) + gamma[zone(a, ond_tsnetwork)] ≤ MOVING_COST * ond_tsnetwork.arccost[a])
    @constraint(model_dual, [a in intersect(nontransfer_moving_arcs,non_f_arc_keys, gamma_arcs)], lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes) + gamma[zone(a, ond_tsnetwork)] <= MOVING_COST * ond_tsnetwork.arccost[a])
    @constraint(model_dual, [a in intersect(nontransfer_moving_arcs,non_f_arc_keys, non_gamma_arcs)], lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes) <= MOVING_COST * ond_tsnetwork.arccost[a])
    @constraint(model_dual, [a in intersect(nontransfer_moving_arcs,f_arc_keys, non_gamma_arcs)], K*beta[a] + lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes)  <= MOVING_COST * ond_tsnetwork.arccost[a])

    @constraint(model_dual, [a in intersect(nontransfer_waiting_arcs, f_arc_keys, gamma_arcs)], K*beta[a] + lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes) + gamma[zone(a, ond_tsnetwork)] <= 0)
    @constraint(model_dual, [a in intersect(nontransfer_waiting_arcs,non_f_arc_keys, gamma_arcs)], lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes) + gamma[zone(a, ond_tsnetwork)] <= 0)
    @constraint(model_dual, [a in intersect(nontransfer_waiting_arcs,non_f_arc_keys, non_gamma_arcs)], lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes) <= 0)
    @constraint(model_dual, [a in intersect(nontransfer_waiting_arcs,f_arc_keys, non_gamma_arcs)], K*beta[a] + lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes)  <= 0)


    return model_dual
end

"""
    retrieve_duals_duals(model, ond_tsnetwork, numcust, non_flow_balance_nodes, moving_arcs, station_waiting_arcs, nontransfer_moving_arcs, nontransfer_waiting_arcs, f_arc_keys, non_f_arc_keys,  gamma_arcs, non_gamma_arcs)
    Once optimized, retrieve the duals of the primal subproblem
"""
function retrieve_duals_duals(model, ond_tsnetwork, numcust, non_flow_balance_nodes, moving_arcs, station_waiting_arcs, nontransfer_moving_arcs, nontransfer_waiting_arcs, f_arc_keys, non_f_arc_keys,  gamma_arcs, non_gamma_arcs)
    if termination_status(model) == MOI.OPTIMAL
        f = Array{Float64}(undef,numcust, ond_tsnetwork.numarcs)
        y = Array{Float64}(undef, ond_tsnetwork.numarcs)
        duals = Float64[]
        for constr in JuMP.all_constraints(model, include_variable_in_set_constraints = true)
            dual_value = JuMP.dual(constr)
            push!(duals, dual_value)
        end

        i= 1
        for n in non_flow_balance_nodes
            i += 1
        end
        for p in 1:numcust
            for a in moving_arcs
                f[p,a] = duals[i]
                i+=1
            end
            for a in station_waiting_arcs
                f[p,a] = duals[i]
                i+=1
            end
        end
        for a in intersect(nontransfer_moving_arcs, f_arc_keys, gamma_arcs)
            y[a] = duals[i]
            i+=1
        end
        for a in intersect(nontransfer_moving_arcs,non_f_arc_keys, gamma_arcs)
            y[a] = duals[i]

            i+=1
        end
        for a in intersect(nontransfer_moving_arcs,non_f_arc_keys, non_gamma_arcs)
            y[a] = duals[i]

            i+=1
        end
        for a in intersect(nontransfer_moving_arcs,f_arc_keys, non_gamma_arcs)
            y[a] = duals[i]

            i+=1
        end

        for a in intersect(nontransfer_waiting_arcs, f_arc_keys, gamma_arcs)
            y[a] = duals[i]

            i+=1
        end
        for a in intersect(nontransfer_waiting_arcs,non_f_arc_keys, gamma_arcs)
            y[a] = duals[i]

            i+=1
        end
        for a in intersect(nontransfer_waiting_arcs,non_f_arc_keys, non_gamma_arcs)
            y[a] = duals[i]

            i+=1
        end
        for a in intersect(nontransfer_waiting_arcs,f_arc_keys, non_gamma_arcs)
            y[a] = duals[i]

            i+=1
        end
    
        return f, y
    else
        println("Le modèle n'a pas été résolu de manière optimale.")
        return 0, 0, 0, 0, 0
    end
end



#Benders algorithm
function Benders(ond_tsnetwork, transit_tsnetwork, numcust, numstations, numzones, horizon, K, W, R, R1, R2, R3, R4, R5, Npo, Npb, Npa, Npd, walking_speed, routeid, routedesc, numroutes, TIME_COSTS, MOVING_COST, WALKING_COST, optimality_gap, time_limit, dual_vs_primal, run_mip_check)
    MP = Model(Gurobi.Optimizer)
    set_optimizer_attributes(
        MP, 
        "TimeLimit" => time_limit, 
        "MIPGap" => optimality_gap, 
        "OutputFlag" => 0)

    #---------------------UTILS----------------------------#
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
    walking_times_to_station = walking_t_to_station(ond_tsnetwork, numstations, walking_speed)
    dist_to_station = distance_to_station(ond_tsnetwork, numstations)

    moving_arcs, waiting_arcs = partition_arcs(ond_tsnetwork)
    flow_balance_nodes = flowBalanceNodes(ond_tsnetwork) #for vehicles

    dist_alight_to_dest = distance_alight_to_dest(ond_tsnetwork, numstations) #Dict[p, s] --> distance from station s to destination of p
    station_waiting_arcs = stationWaitingArcs(ond_tsnetwork, waiting_arcs,numstations)

    d_to_dest = distance_to_dest(ond_tsnetwork)

    f_arc_keys = union(moving_arcs, station_waiting_arcs)
    non_f_arc_keys = []
    for arc in 1:ond_tsnetwork.numarcs
        if !(arc in f_arc_keys)
            push!(non_f_arc_keys, arc)
        end
    end

    nontransfer_arcs = 1:ond_tsnetwork.numarcs

    routeWalkingDict = determineRouteWalkingDistance(ond_tsnetwork, numroutes, routedesc, walking_times_to_station, dist_alight_to_dest)

    boarding_route = boardingNodeRoute(ond_tsnetwork, numroutes, routedesc)
    alighting_route = alightNodeRoute(ond_tsnetwork, numroutes, routedesc)
    reverse_boarding_route = reverseBoardingNodeRoute(boarding_route)
    reverse_alighting_route = reverseAlightNodeRoute(alighting_route)

    r_time_boarding, r_time_alighting = route_boarding_alighting_time(ond_tsnetwork, boarding_route, alighting_route, numroutes)
    r_boarding_station = boardingStationRoute(ond_tsnetwork, numroutes, routedesc)
    r_alighting_station = alightingStationRoute(ond_tsnetwork, numroutes, routedesc)


    

    #------------ SANITY CHECK: THE MIP MODEL -----------------#

    if run_mip_check == true
        f, y , x , zO, zD, preprocessing_t, optimizing_t = model(ond_tsnetwork,numcust, numstations, numzones, horizon,  K, W, R, R1, R2, R3, R4, R5, Npo, Npb, Npa, Npd, WALKING_SPEED, routedesc, numroutes, TIME_COSTS, MOVING_COST, WALKING_COST)
        save("MIP_check.jld", "f", f, "y", y, "x", x, "zO", zO, "zD", zD)
    end


    #--------------------- VARIABLES - Master Problem -------------------------------# 

    @variable(MP, x[p in 1:numcust, r in R[p]], Bin)
    @variable(MP, zO[p in 1:numcust, n in Npo[p]], Bin)
    @variable(MP, zD[p in 1:numcust, n in Npd[p]], Bin)
    @variable(MP, Θ >= 0)


    #--------------------- CONSTRAINTS -------------------------------#

    @constraint(MP, [p in 1:numcust], sum(x[p, r] for r in R[p]) == 1) #26 ; choosing only one route
    @constraint(MP, [p in 1:numcust], sum(zO[p, n] for n in Npo[p]) == sum(x[p,r] for r in R2[p]) + sum(x[p,r] for r in R3[p]) + sum(x[p,r] for r in R5[p])) 
    @constraint(MP, [p in 1:numcust], sum(zD[p, n] for n in Npd[p]) == sum(x[p,r] for r in R2[p]) + sum(x[p,r] for r in R4[p]) + sum(x[p,r] for r in R5[p]))

    @constraint(MP, [p in 1:numcust, n in Npo[p], r in union(R3[p], R5[p])], zO[p,n] * x[p,r] * (r_time_boarding[r] - ond_tsnetwork.nodedesc[n][2] - dist_to_station[p, r_boarding_station[r]]) >= 0)
    @constraint(MP, [p in 1:numcust, n in Npd[p], r in union(R4[p], R5[p])], zD[p,n] * x[p,r] * (ond_tsnetwork.nodedesc[n][2] - r_time_alighting[r] - dist_alight_to_dest[p, r_alighting_station[r]]) >= 0)

    @constraint(MP, [p in 1:numcust, no in Npo[p], nd in Npd[p]],  zO[p, no] * zD[p, nd] * (ond_tsnetwork.nodedesc[nd][2] - ond_tsnetwork.nodedesc[no][2] - d_to_dest[p])  >= 0)


    @objective(MP, Min,
    sum(TIME_COSTS[ond_tsnetwork.nodedesc[node][2]]*zO[p, node] for p in 1:numcust, node in Npo[p]) 
    + sum(TIME_COSTS[ond_tsnetwork.nodedesc[node][2]]*zD[p, node] for p in 1:numcust, node in Npd[p])
    + sum(TIME_COSTS[arrivalTimesType13(ond_tsnetwork, numroutes, numstations, r, routedesc, walking_speed, horizon, p)] * x[p,r] for p in 1:numcust, r in union(R1[p], R3[p]))
    + WALKING_COST * sum( routeWalkingDict[p, r] * x[p,r] for p in 1:numcust, r in R[p])
    + Θ
    )  

    lower_bound_all = []
    upper_bound_all = []
    MP_time = []
    SP_time = []
    it = 1
    while true
        #Master Problem
        println("ITERATION: ", it)
        push!(MP_time, @elapsed optimize!(MP))
        lower_bound_new = objective_value(MP)
        push!(lower_bound_all, lower_bound_new)
        x_val, zO_val, zD_val, Θ_val = value.(MP[:x]), value.(MP[:zO]), value.(MP[:zD]), value.(MP[:Θ])
        save("Benders/master_problem_$it.jld", "x", x_val, "zO", zO_val, "zD", zD_val, "Θ", Θ_val)
        println("Master problem solution value ", lower_bound_new)


        #Dual Subproblem

        if dual_vs_primal == "DUAL" #In this case we solve directly the dual
            SP_dual = subproblem_dual(ond_tsnetwork, numcust, numzones, K, MOVING_COST, optimality_gap, moving_arcs, waiting_arcs, nontransfer_arcs, flow_balance_nodes, f_arc_keys, non_f_arc_keys, station_waiting_arcs)
            push!(SP_time, @elapsed optimize!(SP_dual))
            println("Dual termination status : ", termination_status(SP_dual))
            α, β, λ, γ, ξ = value.(SP_dual[:alpha]), value.(SP_dual[:beta]), value.(SP_dual[:lambda]), value.(SP_dual[:gamma]), value.(SP_dual[:ksi])

        
        elseif dual_vs_primal == "PRIMAL" #In this case we solve the primal and retrieve the duals
            SP_dual = subproblem_primal(ond_tsnetwork, numcust, numstations, K, W, R3, R4, R5, Npo, Npb, Npa, Npd, WALKING_SPEED, routedesc, numroutes, MOVING_COST, x_val, zO_val, zD_val, time_limit)
            println("check1")
            push!(SP_time, @elapsed optimize!(SP_dual))
            println("SP termination status : ", termination_status(SP_dual))
            α, β, λ, γ, ξ = retrieve_duals(SP_dual, ond_tsnetwork, numcust, Npo, Npb, Npa, Npd, other_nodes, flow_balance_nodes, nontransfer_arcs, f_arc_keys, I)

        end

        #Infeasibility cuts
        if termination_status(SP_dual) == MOI.INFEASIBLE
            @constraint(
                MP, 
                0 
                ≥   sum( sum(α[p, n] * zO[p,n] for n in Npo[p])
                - sum(α[p,n] * sum(x[p,r] for r in intersect(union(R3[p], R5[p]), reverse_boarding_route[n])) for n in Npb[p])
                + sum(α[p,n] *  sum(x[p,r] for r in intersect(union(R4[p], R5[p]), reverse_alighting_route[n]))  for n in Npa[p])
                - sum(α[p, n] * zD[p,n] for n in Npd[p]) for p in 1:numcust) 
            + sum(W[k] * γ[k] for k in 1:numzones) 
            - sum(ξ[p,a] for p in 1:numcust for a in f_arc_keys)
            )
            global obj_SP = Inf
            println("SP INFEASIBLE. Adding feasibility cut...")

        # If dual subproblem is bounded and solves to optimality, add optimality cut
        elseif termination_status(SP_dual) == MOI.OPTIMAL
            global obj_SP = objective_value(SP_dual)
            #SP_primal = subproblem_primal(ond_tsnetwork, numcust, numstations, K, W, R3, R4, R5, Npo, Npb, Npa, Npd, WALKING_SPEED, routedesc, numroutes, MOVING_COST, x, zO, zD, time_limit)
            #println("CHECK | Primal value = ", primal_value)     
            println("Adding optimality cut...")
            @constraint(MP, Θ >=  sum( sum(α[p, n] * zO[p,n] for n in Npo[p])
            - sum(α[p,n] * sum(x[p,r] for r in intersect(union(R3[p], R5[p]), reverse_boarding_route[n])) for n in Npb[p])
            + sum(α[p,n] *  sum(x[p,r] for r in intersect(union(R4[p], R5[p]), reverse_alighting_route[n]))  for n in Npa[p])
            - sum(α[p, n] * zD[p,n] for n in Npd[p]) for p in 1:numcust) 
        + sum(W[k] * γ[k] for k in 1:numzones) 
        - sum(ξ[p,a] for p in 1:numcust for a in f_arc_keys))
        end

        println("Subproblem solution value ", obj_SP)
        upper_bound_new = lower_bound_new + obj_SP - Θ_val
        push!(upper_bound_all,  upper_bound_new)

        # Termination criteria

        if (sum(MP_time) + sum(SP_time) ≥ time_limit)
            println("Time limit reached")
            break
        elseif (abs( (upper_bound_new - lower_bound_new) / lower_bound_new) < optimality_gap)
            println("Optimality gap reached")
            println(abs(upper_bound_new - lower_bound_new) / lower_bound_new)
            break
        end
        it += 1

        println("------------------------------------------------------------------")
        println("Upper bounds: ", upper_bound_all)
        println("Lower bounds: ", lower_bound_all)
        println("MP times (s): ", MP_time)
        println("SP times (s): ", SP_time)
        println("------------------------------------------------------------------")

    end
    return upper_bound_all, lower_bound_all, MP_time, SP_time

end

function Benders_pareto(ond_tsnetwork, transit_tsnetwork, numcust, numstations, numzones, horizon, K, W, R, R1, R2, R3, R4, R5, Npo, Npb, Npa, Npd, walking_speed, routeid, routedesc, numroutes, TIME_COSTS, MOVING_COST, WALKING_COST, optimality_gap, time_limit, dual_vs_primal, run_mip_check)
    MP = Model(Gurobi.Optimizer)
    set_optimizer_attributes(
        MP, 
        "TimeLimit" => time_limit, 
        "MIPGap" => optimality_gap, 
        "OutputFlag" => 0)

    #---------------------UTILS----------------------------#
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
    walking_times_to_station = walking_t_to_station(ond_tsnetwork, numstations, walking_speed)
    dist_to_station = distance_to_station(ond_tsnetwork, numstations)

    moving_arcs, waiting_arcs = partition_arcs(ond_tsnetwork)
    flow_balance_nodes = flowBalanceNodes(ond_tsnetwork) #for vehicles

    dist_alight_to_dest = distance_alight_to_dest(ond_tsnetwork, numstations) #Dict[p, s] --> distance from station s to destination of p
    station_waiting_arcs = stationWaitingArcs(ond_tsnetwork, waiting_arcs,numstations)

    d_to_dest = distance_to_dest(ond_tsnetwork)

    f_arc_keys = union(moving_arcs, station_waiting_arcs)
    non_f_arc_keys = []
    for arc in 1:ond_tsnetwork.numarcs
        if !(arc in f_arc_keys)
            push!(non_f_arc_keys, arc)
        end
    end

    nontransfer_arcs = 1:ond_tsnetwork.numarcs

    routeWalkingDict = determineRouteWalkingDistance(ond_tsnetwork, numroutes, routedesc, walking_times_to_station, dist_alight_to_dest)

    boarding_route = boardingNodeRoute(ond_tsnetwork, numroutes, routedesc)
    alighting_route = alightNodeRoute(ond_tsnetwork, numroutes, routedesc)
    reverse_boarding_route = reverseBoardingNodeRoute(boarding_route)
    reverse_alighting_route = reverseAlightNodeRoute(alighting_route)

    r_time_boarding, r_time_alighting = route_boarding_alighting_time(ond_tsnetwork, boarding_route, alighting_route, numroutes)
    r_boarding_station = boardingStationRoute(ond_tsnetwork, numroutes, routedesc)
    r_alighting_station = alightingStationRoute(ond_tsnetwork, numroutes, routedesc)

    t1s = arc_origin_times(ond_tsnetwork)
    gamma_arcs = Set()
    for arc in 1:ond_tsnetwork.numarcs
        if zone(arc, ond_tsnetwork) in 1:numzones 
            if t1s[arc] == 0
                push!(gamma_arcs, arc)
            end
        end
    end
    non_gamma_arcs = []
    for arc in 1:ond_tsnetwork.numarcs
        if !(arc in gamma_arcs)
            push!(non_gamma_arcs, arc)
        end
    end
    nontransfer_moving_arcs = intersect(moving_arcs, nontransfer_arcs)
    nontransfer_waiting_arcs = intersect(waiting_arcs, nontransfer_arcs)

    non_flow_balance_nodes = []
    for n in 1:ond_tsnetwork.numnodes
        if !(n in flow_balance_nodes)
            push!(non_flow_balance_nodes, n)
        end
    end


    #------------ SANITY CHECK: THE MIP MODEL -----------------#

    if run_mip_check == true
        f, y , x , zO, zD, preprocessing_t, optimizing_t = model(ond_tsnetwork,numcust, numstations, numzones, horizon,  K, W, R, R1, R2, R3, R4, R5, Npo, Npb, Npa, Npd, WALKING_SPEED, routedesc, numroutes, TIME_COSTS, MOVING_COST, WALKING_COST)
        save("MIP_check.jld", "f", f, "y", y, "x", x, "zO", zO, "zD", zD)
    end


    #--------------------- VARIABLES -------------------------------# 

    @variable(MP, x[p in 1:numcust, r in R[p]], Bin)
    @variable(MP, zO[p in 1:numcust, n in Npo[p]], Bin)
    @variable(MP, zD[p in 1:numcust, n in Npd[p]], Bin)
    @variable(MP, Θ >= 0)


    #--------------------- CONSTRAINTS -------------------------------#

    @constraint(MP, [p in 1:numcust], sum(x[p, r] for r in R[p]) == 1) #26 ; choosing only one route
    @constraint(MP, [p in 1:numcust], sum(zO[p, n] for n in Npo[p]) == sum(x[p,r] for r in R2[p]) + sum(x[p,r] for r in R3[p]) + sum(x[p,r] for r in R5[p])) 
    @constraint(MP, [p in 1:numcust], sum(zD[p, n] for n in Npd[p]) == sum(x[p,r] for r in R2[p]) + sum(x[p,r] for r in R4[p]) + sum(x[p,r] for r in R5[p]))

    @constraint(MP, [p in 1:numcust, n in Npo[p], r in union(R3[p], R5[p])], zO[p,n] * x[p,r] * (r_time_boarding[r] - ond_tsnetwork.nodedesc[n][2] - dist_to_station[p, r_boarding_station[r]]) >= 0)
    @constraint(MP, [p in 1:numcust, n in Npd[p], r in union(R4[p], R5[p])], zD[p,n] * x[p,r] * (ond_tsnetwork.nodedesc[n][2] - r_time_alighting[r] - dist_alight_to_dest[p, r_alighting_station[r]]) >= 0)

    @constraint(MP, [p in 1:numcust, no in Npo[p], nd in Npd[p]],  zO[p, no] * zD[p, nd] * (ond_tsnetwork.nodedesc[nd][2] - ond_tsnetwork.nodedesc[no][2] - d_to_dest[p])  >= 0)


    @objective(MP, Min,
    sum(TIME_COSTS[ond_tsnetwork.nodedesc[node][2]]*zO[p, node] for p in 1:numcust, node in Npo[p]) 
    + sum(TIME_COSTS[ond_tsnetwork.nodedesc[node][2]]*zD[p, node] for p in 1:numcust, node in Npd[p])
    + sum(TIME_COSTS[arrivalTimesType13(ond_tsnetwork, numroutes, numstations, r, routedesc, walking_speed, horizon, p)] * x[p,r] for p in 1:numcust, r in union(R1[p], R3[p]))
    + WALKING_COST * sum( routeWalkingDict[p, r] * x[p,r] for p in 1:numcust, r in R[p])
    + Θ
    )  

    lower_bound_all = []
    upper_bound_all = []
    MP_time = []
    SP_time = []
    Pareto_time = []
    it = 1
    while true 
        #Master Problem
        println("ITERATION: ", it)
        push!(MP_time, @elapsed optimize!(MP))
        lower_bound_new = objective_value(MP)
        push!(lower_bound_all, lower_bound_new)
        x_val, zO_val, zD_val, Θ_val = value.(MP[:x]), value.(MP[:zO]), value.(MP[:zD]), value.(MP[:Θ])
        save("Benders/pareto_master_problem_$it.jld", "x", x_val, "zO", zO_val, "zD", zD_val, "Θ", Θ_val)
        println("Master problem solution value ", lower_bound_new)


        #Dual Subproblem

        if dual_vs_primal == "DUAL"
            SP_dual = subproblem_dual(ond_tsnetwork, numcust, numzones, K, MOVING_COST, optimality_gap, moving_arcs, waiting_arcs, nontransfer_arcs, flow_balance_nodes, f_arc_keys, non_f_arc_keys, station_waiting_arcs)
            push!(SP_time, @elapsed optimize!(SP_dual))
            α, β, λ, γ, ξ = value.(SP_dual[:alpha]), value.(SP_dual[:beta]), value.(SP_dual[:lambda]), value.(SP_dual[:gamma]), value.(SP_dual[:ksi])
            println("Dual termination status : ", termination_status(SP_dual))
        
        elseif dual_vs_primal == "PRIMAL" 
            SP_dual = subproblem_primal(ond_tsnetwork, numcust, numstations, K, W, R3, R4, R5, Npo, Npb, Npa, Npd, WALKING_SPEED, routedesc, numroutes, MOVING_COST, x_val, zO_val, zD_val, time_limit)
            push!(SP_time, @elapsed optimize!(SP_dual))
            println("SP termination status : ", termination_status(SP_dual))
        end

        if termination_status(SP_dual) == MOI.INFEASIBLE
            α, β, λ, γ, ξ = retrieve_duals(SP_dual, ond_tsnetwork, numcust, Npo, Npb, Npa, Npd, other_nodes, flow_balance_nodes, nontransfer_arcs, f_arc_keys, I)
            println("SP INFEASIBLE. Adding feasibility cut...")
            @constraint(
                MP, 
                0 
                ≥   sum( sum(α[p, n] * zO[p,n] for n in Npo[p])
                - sum(α[p,n] * sum(x[p,r] for r in intersect(union(R3[p], R5[p]), reverse_boarding_route[n])) for n in Npb[p])
                + sum(α[p,n] *  sum(x[p,r] for r in intersect(union(R4[p], R5[p]), reverse_alighting_route[n]))  for n in Npa[p])
                - sum(α[p, n] * zD[p,n] for n in Npd[p]) for p in 1:numcust) 
            + sum(W[k] * γ[k] for k in 1:numzones) 
            - sum(ξ[p,a] for p in 1:numcust for a in f_arc_keys)
            )
            global obj_SP = Inf
            push!(Pareto_time, 0)


        # If dual subproblem is bounded and solves to optimality, add optimality cut after having solved the Pareto model
        elseif termination_status(SP_dual) == MOI.OPTIMAL
            global obj_SP = objective_value(SP_dual)

            paretoModel = Model(Gurobi.Optimizer)
            set_optimizer_attributes(
                paretoModel, 
                "OutputFlag" => 1, 
                "InfUnbdInfo" => 1)
        

            #--------------------- VARIABLES - Pareto Model -------------------------------# 

            @variable(paretoModel, -1000 <= alpha[p in 1:numcust, n in 1:ond_tsnetwork.numnodes] <= 10)
            @variable(paretoModel, 0<= beta[a in f_arc_keys])
            @variable(paretoModel, lambda[n in 1:ond_tsnetwork.numnodes])
            @variable(paretoModel, gamma[k in 1:numzones])
            @variable(paretoModel, 0<= ksi[p in 1:numcust, a in f_arc_keys])

            #--------------------- CONSTRAINTS -------------------------------#
            @constraint(paretoModel, [n in non_flow_balance_nodes], lambda[n] == 0)

            @constraint(paretoModel, [p in 1:numcust, a in moving_arcs], alpha[p, ond_tsnetwork.arcdesc[a][1]] - alpha[p, ond_tsnetwork.arcdesc[a][2]] - beta[a] - ksi[p, a] <= MOVING_COST * ond_tsnetwork.arccost[a])
            @constraint(paretoModel, [p in 1:numcust, a in station_waiting_arcs], alpha[p, ond_tsnetwork.arcdesc[a][1]] - alpha[p, ond_tsnetwork.arcdesc[a][2]] - beta[a] - ksi[p, a] <= 0)

            @constraint(paretoModel, [a in intersect(nontransfer_moving_arcs, f_arc_keys, gamma_arcs)], K*beta[a] + lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes) + gamma[zone(a, ond_tsnetwork)] ≤ MOVING_COST * ond_tsnetwork.arccost[a])
            @constraint(paretoModel, [a in intersect(nontransfer_moving_arcs,non_f_arc_keys, gamma_arcs)], lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes) + gamma[zone(a, ond_tsnetwork)] <= MOVING_COST * ond_tsnetwork.arccost[a])
            @constraint(paretoModel, [a in intersect(nontransfer_moving_arcs,non_f_arc_keys, non_gamma_arcs)], lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes) <= MOVING_COST * ond_tsnetwork.arccost[a])
            @constraint(paretoModel, [a in intersect(nontransfer_moving_arcs,f_arc_keys, non_gamma_arcs)], K*beta[a] + lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes)  <= MOVING_COST * ond_tsnetwork.arccost[a])

            @constraint(paretoModel, [a in intersect(nontransfer_waiting_arcs, f_arc_keys, gamma_arcs)], K*beta[a] + lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes) + gamma[zone(a, ond_tsnetwork)] <= 0)
            @constraint(paretoModel, [a in intersect(nontransfer_waiting_arcs,non_f_arc_keys, gamma_arcs)], lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes) + gamma[zone(a, ond_tsnetwork)] <= 0)
            @constraint(paretoModel, [a in intersect(nontransfer_waiting_arcs,non_f_arc_keys, non_gamma_arcs)], lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes) <= 0)
            @constraint(paretoModel, [a in intersect(nontransfer_waiting_arcs,f_arc_keys, non_gamma_arcs)], K*beta[a] + lambda[ond_tsnetwork.arcdesc[a][1]]*indicator(ond_tsnetwork.arcdesc[a][1], flow_balance_nodes) - lambda[ond_tsnetwork.arcdesc[a][2]]*indicator(ond_tsnetwork.arcdesc[a][2], flow_balance_nodes)  <= 0)
        

            @constraint(paretoModel, c_pareto, 
            obj_SP - 0.00001  <= sum(sum(alpha[p, n] * zO_val[p,n] for n in Npo[p])
            - sum(alpha[p,n] * sum(x_val[p,r] for r in intersect(union(R3[p], R5[p]), reverse_boarding_route[n])) for n in Npb[p])
            + sum(alpha[p,n] *  sum(x_val[p,r] for r in intersect(union(R4[p], R5[p]), reverse_alighting_route[n]))  for n in Npa[p])
            - sum(alpha[p, n] * zD_val[p,n] for n in Npd[p]) for p in 1:numcust) 
            + sum(W[k] * gamma[k] for k in 1:numzones) - sum(ksi[p,a] for p in 1:numcust, a in f_arc_keys) <= obj_SP)

            @objective(paretoModel, Max, sum(sum(alpha[p, n] * 0.5 for n in Npo[p])
            - sum(alpha[p,n] * sum(0.5 for r in intersect(union(R3[p], R5[p]), reverse_boarding_route[n])) for n in Npb[p])
            + sum(alpha[p,n] *  sum(0.5 for r in intersect(union(R4[p], R5[p]), reverse_alighting_route[n]))  for n in Npa[p])
            - sum(alpha[p, n] * 0.5 for n in Npd[p]) for p in 1:numcust) 
            + sum(W[k] * gamma[k] for k in 1:numzones) - sum(ksi[p,a] for p in 1:numcust, a in f_arc_keys))

            println("Compiling Pareto optimal cut...")
            push!(Pareto_time, @elapsed optimize!(paretoModel))
            println("Pareto value ", objective_value(paretoModel))
            println("Pareto status ", termination_status(paretoModel))

            α, β, λ, γ, ξ = value.(paretoModel[:alpha]), value.(paretoModel[:beta]), value.(paretoModel[:lambda]), value.(paretoModel[:gamma]), value.(paretoModel[:ksi])

            println("Adding optimality cut...")
            @constraint(MP, Θ >=  sum( sum(α[p, n] * zO[p,n] for n in Npo[p])
            - sum(α[p,n] * sum(x[p,r] for r in intersect(union(R3[p], R5[p]), reverse_boarding_route[n])) for n in Npb[p])
            + sum(α[p,n] *  sum(x[p,r] for r in intersect(union(R4[p], R5[p]), reverse_alighting_route[n]))  for n in Npa[p])
            - sum(α[p, n] * zD[p,n] for n in Npd[p]) for p in 1:numcust) 
            + sum(W[k] * γ[k] for k in 1:numzones) 
            - sum(ξ[p,a] for p in 1:numcust for a in f_arc_keys))

        end

        println("Subproblem solution value ", obj_SP)
        upper_bound_new = lower_bound_new + obj_SP - Θ_val
        push!(upper_bound_all,  upper_bound_new)

        # Termination criteria

        if (sum(MP_time) + sum(SP_time) + sum(Pareto_time) ≥ time_limit)
            println("Time limit reached")
            break
        elseif (abs( (upper_bound_new - lower_bound_new) / lower_bound_new) < optimality_gap)
            println("Optimality gap reached")
            println(abs(upper_bound_new - lower_bound_new) / lower_bound_new)
            break
        end
        it += 1

        println("------------------------------------------------------------------")
        println("Upper bounds: ", upper_bound_all)
        println("Lower bounds: ", lower_bound_all)
        println("MP times (s): ", MP_time)
        println("SP times (s): ", SP_time)
        println("Pareto times (s): ", Pareto_time)
        println("------------------------------------------------------------------")

    end
    return upper_bound_all, lower_bound_all, MP_time, SP_time, Pareto_time

end



function runBenders(optimality_gap, time_limit, primal_vs_dual, pareto, run_mip_check)
    setting = "CROSS" #"CROSS" or "LINEAR_EXPANSION"

#----------------------------------PARAMETERS----------------------------------#  	

    SAMPLE_SIZE = 3
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


    DEPARTURE_INTERVAL = 5         #Size of the departure window
    departure_times = rand(1:DEPARTURE_INTERVAL, numcust)

    DELTA = 20
    DELTA2 = DELTA
    iter = 1
    RADIUS = 5

    if setting == "CROSS"
        run(`/Users/Meilame/opt/anaconda3/bin/python generate_cross.py $RADIUS $SAMPLE_SIZE $DELTA $DELTA2 $iter`)
    elseif setting == "LINEAR_EXPANSION"
        run(`/Users/Meilame/opt/anaconda3/bin/python generate_linear.py $RADIUS $SAMPLE_SIZE $DELTA $DELTA2 $iter`)
    end

    transit_tsnetwork = transit_createfullnetwork(transit_locationfilename, transit_arcfilename, numstations, horizon, tstep, freq, MAX_WAITING_TRANSFER, setting)
    ond_tsnetwork = ond_createfullnetwork(ond_locationfilename, ond_arcfilename, transit_tsnetwork, numlocs, numcust, numzones, horizon, tstep)

    R, R1, R2, R3, R4, R5, Npo, Npb, Npa, Npd, routeid, routedesc, numroutes, walking_dict = preprocessing(ond_tsnetwork, transit_tsnetwork, numstations, MAX_WALKING, MAX_WAITING, WALKING_SPEED, departure_times)
    
    if pareto == false
        upper_bound_all, lower_bound_all, MP_time, SP_time = Benders(ond_tsnetwork, transit_tsnetwork, numcust, numstations, numzones, horizon, K, W, R, R1, R2, R3, R4, R5, Npo, Npb, Npa, Npd, WALKING_SPEED, routeid, routedesc, numroutes, TIME_COSTS, MOVING_COST, WALKING_COST, optimality_gap, time_limit, primal_vs_dual, run_mip_check)
    else 
        upper_bound_all, lower_bound_all, MP_time, SP_time, Pareto_time = Benders_pareto(ond_tsnetwork, transit_tsnetwork, numcust, numstations, numzones, horizon, K, W, R, R1, R2, R3, R4, R5, Npo, Npb, Npa, Npd, WALKING_SPEED, routeid, routedesc, numroutes, TIME_COSTS, MOVING_COST, WALKING_COST, optimality_gap, time_limit, primal_vs_dual, run_mip_check)
    end

    println("Upper bound = ", upper_bound_all)
    println("Lower bound = ", lower_bound_all)


    plotting_UB = []
    push!(plotting_UB, upper_bound_all[1])
    for i in 2:length(upper_bound_all)
        push!(plotting_UB, min(upper_bound_all[i], plotting_UB[i-1]))
    end
    
    plot(
    [plotting_UB lower_bound_all], 
    label = ["Upper bound" "Lower bound"], 
    xlabel = "Iteration",
    ylabel = "Value")

    savefig("iterations.png")

    plot(
    cumsum(MP_time .+ SP_time), 
    [plotting_UB lower_bound_all], 
    label = ["Upper bound" "Lower bound"], 
    xlabel = "Time (s)",
    ylabel = "Value",
)

    savefig("time.png")
end

#------------ RUN BENDERS ALGORITHM ------------#
#Comment the last line of run_optimization.jl before running !

optimality_gap = 1e-4 #Increase it makes the algorithm converge faster (but solution might be far from optimum), reducing it makes the solution better (but increases runtime)
time_limit = 3600
pareto = true #true to use Pareto-optimal cuts, false not to
run_mip_check = false #First run the MILP model to perform sanity check at the end? Only enable for small size settings ! otherwise the MILP itself will take long

runBenders(optimality_gap, time_limit, "PRIMAL", pareto, run_mip_check)
