include("energy_system.jl")
include("input_data.jl")

# This script contains models
################### Models ####################

################### PDE MODEL ####################
function build_PDE_model(ES::GasSystem)

    nodes, sources, demands, pipes, compressors,
        N, S, D, P, T, C, dt, C_cur_load_gas = get_gas_system_attributes(ES)

    ut1, ut2 = PDE_type_parameters(ES.PDE_type)

    #################### Initializing model ####################
    m = Model()
    m = set_solver(m, ES)

    #################### Variables ####################
    @variables m begin
        Q_w[t in T, s in S] # Scheduled gas source flowrate [kg/s] during each timestep
        pr[t in T, n in N] # Nodal pressure  [MPa]
        Q_nm[t in T, p in P] # Gas flow in each pipe [kg/s] during each timestep
        Q_nm_in[t in T, p in P] # Gas flow entering the pipe
        Q_nm_out[t in T, p in P] # Gas flow exiting the pipe
        pr_avg[t in T, p in P] # Average pressure in each pipe (MPa)
        Qd_cur_gas[t in T, d in D] # load shedding at node n at timestep t [kg/s]
    end
    
    if isempty(compressors) == false # if there are compressors in the nework
        @variable(m, Q_c[t in T, c in C]) # Gas flow in each compressor [kg/s]
    end

    ################### Objective function #################### Why is it quadratic?
    if ES.obj_type == "linear"
        @objective(m, Min, dt/3600 *(
            sum(Q_w[t,s]*sources[s].c1 for s in S, t in T) +
            sum(Qd_cur_gas[t,d]*C_cur_load_gas for d in D, t in T) +
            (isempty(compressors) ? 0.0 : sum(compressors[c].compr_cost*(
                pr[t,compressors[c].stop]-pr[t,compressors[c].start]) for c in C, t in T))
            ))
    elseif ES.obj_type == "quadratic"
        @objective(m, Min, dt/3600 *(
            sum(Q_w[t,s]*sources[s].c1 + Q_w[t,s]^2*sources[s].c2 for s in S, t in T) +
            sum(Qd_cur_gas[t,d]*C_cur_load_gas for d in D, t in T) +
            (isempty(compressors) ? 0.0 : sum(compressors[c].compr_cost*(
                pr[t,compressors[c].stop]-pr[t,compressors[c].start]) for c in C, t in T))
            ))
    end

    ################### Constraints ####################
    @constraints m begin
        supply_limits[t in T, s in S],
            sources[s].Q_min <= Q_w[t,s] <= sources[s].Q_max
        # Only the pressure for the original nodes must be constrained
        pressure_limits[t in T, n in N],
            nodes[n].pr_min <= pr[t,n] <= nodes[n].pr_max
        curtailment_limit[t in T, d in D],
            0 <= Qd_cur_gas[t,d] <= demands[d].Qd_gas[t]
        linepack_end_condition, # Restore linepack at the end of the time horizon
            sum(Q_w[t,s] for t in T, s in S) >= sum(demands[d].Qd_gas[t] for t in T, d in D)
    end

    @expression(m, expr_node_balance[t in T, n in N],
        sum(Q_w[t,s] for s in S if sources[s].location == n) -
        sum(demands[d].Qd_gas[t] for d in D if demands[d].node_load == n) +
        sum(Qd_cur_gas[t,d] for d in D if demands[d].node_load == n) +
        sum(Q_nm_out[t,p] for p in P if pipes[p].stop == n) -
        sum(Q_nm_in[t,p] for p in P if pipes[p].start == n))

    # Nodal balance equation with and without compressors
    if isempty(compressors) == true # if there are no compressors in the nework
        @constraint(m, node_balance[t in T, n in N], expr_node_balance[t,n] == 0)
    else 
        @constraints m begin
            node_balance[t in T, n in N],
                expr_node_balance[t,n] +
                sum(Q_c[t,c] for c in C if compressors[c].stop == n) -
                sum(Q_c[t,c] for c in C if compressors[c].start == n ) == 0
            compressor_flow[t in T, c in C],
                0 <= Q_c[t,c]
            compressor_limits_l[t in T, c in C],
                compressors[c].CR_min*pr[t,compressors[c].start]  <= pr[t,compressors[c].stop] 
            compressor_limits_u[t in T, c in C],
                pr[t,compressors[c].stop] <= compressors[c].CR_max*pr[t,compressors[c].start]
        end
    end

    # Gas flow model nonconvex
    @constraints m begin
        avg_flow_def[t in T, p in P], 
            Q_nm[t,p] == (Q_nm_in[t,p]+Q_nm_out[t,p])/2
        avg_pressure_def[t in T, p in P],
            pr_avg[t,p] == (pr[t,pipes[p].start]+pr[t,pipes[p].stop])/2
        conservation_mass_initial[t=[1], p in P],
            pipes[p].V_m*(Q_nm_out[t,p]-Q_nm_in[t,p])/pipes[p].LL == 0
        conservation_mass[t in T, p in P; t >= 2],
            ut1*(pr_avg[t,p]-pr_avg[t-1,p])/dt +
            pipes[p].V_m*(Q_nm_out[t,p]-Q_nm_in[t,p])/pipes[p].LL == 0
    end

    @NLconstraints m begin
        conservation_momentum_initial[t=[1], p in P],
            pr_avg[t,p]*pipes[p].V_p*(pr[t,pipes[p].stop]-pr[t,pipes[p].start])/pipes[p].LL +
            pipes[p].V_f*(abs(Q_nm[t,p])*Q_nm[t,p]) == 0
        conservation_momentum[t in T, p in P; t >= 2],
            pr_avg[t,p]*ut2*(Q_nm[t,p]-Q_nm[t-1,p])/dt +
            pr_avg[t,p]*pipes[p].V_p*(pr[t,pipes[p].stop]-pr[t,pipes[p].start])/pipes[p].LL +
            pipes[p].V_f*(abs(Q_nm[t,p])*Q_nm[t,p]) == 0
    end
    
    if ES.print_model == 1
        print(m)
    end

    return m
end


################### PDE MODEL ####################
function build_PDE_model_SLP(ES::GasSystem)

    nodes, sources, demands, pipes, compressors,
        N, S, D, P, T, C, dt, C_cur_load_gas = get_gas_system_attributes(ES)

    ut1, ut2 = PDE_type_parameters(ES.PDE_type)

    #################### Initializing model ####################
    m = Model()
    m = set_solver(m, ES)

    #################### Variables ####################
    @variables m begin
        Q_w[t in T, s in S] # Scheduled gas source flowrate [kg/s] during each timestep
        pr[t in T, n in N] # Nodal pressure  [MPa]
        pr_avg[t in T, p in P] # Average pressure in each pipe (MPa)
        Q_nm[t in T, p in P] # Gas flow in each pipe [kg/s] during each timestep
        Q_nm_in[t in T, p in P] # Gas flow entering the pipe
        Q_nm_out[t in T, p in P] # Gas flow exiting the pipe
        Qd_cur_gas[t in T, d in D] # load shedding at node n at timestep t [kg/s]
        γ[t in T, p in P] # Auxiliary variable Taylor Series Expansion
    end

    if isempty(compressors)==false # if there are compressors in the nework
        @variable(m, Q_c[t in T, c in C]) # Gas flow in each compressor [kg/s]
    end

    ################### Objective function ####################
    # Is intitialized within algorithm

    ################### Constraints ####################
    @constraints m begin
        supply_limits[t in T, s in S],
            sources[s].Q_min <= Q_w[t,s] <= sources[s].Q_max
        pressure_limits[t in T, n in N],
            nodes[n].pr_min <= pr[t,n] <= nodes[n].pr_max
        curtailment_limit[t in T, d in D],
            0 <= Qd_cur_gas[t,d] <= demands[d].Qd_gas[t]
        linepack_end_condition, # Restore linepack at the end of the time horizon
            sum(Q_w[t,s] for t in T, s in S) >= sum(demands[d].Qd_gas[t] for t in T, d in D)
    end

    @expression(m, expr_node_balance[t in T, n in N],
        sum(Q_w[t,s] for s in S if sources[s].location == n) -
        sum(demands[d].Qd_gas[t] for d in D if demands[d].node_load == n) +
        sum(Qd_cur_gas[t,d] for d in D if demands[d].node_load == n) +
        sum(Q_nm_out[t,p] for p in P if pipes[p].stop == n) -
        sum(Q_nm_in[t,p] for p in P if pipes[p].start == n))

    # Nodal balance equation with and without compressors
    if isempty(compressors) == true # if there are no compressors in the nework
        @constraint(m, node_balance[t in T, n in N], expr_node_balance[t,n] == 0)
    else 
        @constraints m begin
            node_balance[t in T, n in N],
                expr_node_balance[t,n] +
                sum(Q_c[t,c] for c in C if compressors[c].stop == n) -
                sum(Q_c[t,c] for c in C if compressors[c].start == n ) == 0
            compressor_flow[t in T, c in C],
                0 <= Q_c[t,c]
            compressor_limits_l[t in T, c in C],
                compressors[c].CR_min*pr[t,compressors[c].start] <= pr[t,compressors[c].stop] 
            compressor_limits_u[t in T, c in C],
                pr[t,compressors[c].stop] <= compressors[c].CR_max*pr[t,compressors[c].start]
        end
    end

    # Gas flow model nonconvex
    @constraints m begin
        avg_flow_def[t in T, p in P],
            Q_nm[t,p] == (Q_nm_in[t,p]+Q_nm_out[t,p])/2
        avg_pressure_def[t in T, p in P],
            pr_avg[t,p] == (pr[t,pipes[p].start]+pr[t,pipes[p].stop])/2
        conservation_mass_initial[t=[1], p in P],
            pipes[p].V_m*(Q_nm_out[t,p]-Q_nm_in[t,p])/pipes[p].LL == 0
        conservation_mass[t in T, p in P; t >= 2],
            ut1*(pr_avg[t,p]-pr_avg[t-1,p])/dt +
            pipes[p].V_m*(Q_nm_out[t,p]-Q_nm_in[t,p])/pipes[p].LL == 0
    end

    @constraints m begin
        conservation_momentum_initial[t=[1], p in P],
            pipes[p].V_p * (pr[t,pipes[p].stop]-pr[t,pipes[p].start])/pipes[p].LL +
            pipes[p].V_f * γ[t,p] == 0
        conservation_momentum[t in T, p in P; t >= 2],
            ut2 * (Q_nm[t,p]-Q_nm[t-1,p])/dt +
            pipes[p].V_p * (pr[t,pipes[p].stop]-pr[t,pipes[p].start])/pipes[p].LL +
            pipes[p].V_f * γ[t,p] == 0
    end

    if ES.print_model == 1
        print(m)
    end

    return m
end


# ####################### MISOCP #########################
function build_MISOCP_model(ES)

    nodes, sources, demands, pipes, compressors,
        N, S, D, P, T, C, dt, C_cur_load_gas = get_gas_system_attributes(ES)
    M = ES.M

    m = Model()
    m = set_solver(m, ES)
 
    #################### Variables ####################
    @variables m begin
        Q_w[t in T, s in S] # Scheduled gas source flowrate [kg] during each timestep
        pr[t in T, n in N] # Nodal pressure  [MPa]
        ϕ⁺[t in T, p in P] # auxiliary variables for pressures ϕ⁺ = pr_m + pr_n
        ϕ⁻[t in T, p in P] # auxiliary variables for pressures ϕ⁻ = pr_m - pr_n
        Q_nm[t in T, p in P] # bidirectional gas flow
        Q_nm⁺[t in T, p in P] >= 0 # Gas flow in each pipe [kg] during each timestep
        Q_nm_in⁺[t in T, p in P] >= 0 # Gas flow entering the pipe
        Q_nm_out⁺[t in T, p in P] >= 0 # Gas flow exiting the pipe
        Q_nm⁻[t in T, p in P] >= 0 # Gas flow in each pipe [kg] during each timestep
        Q_nm_in⁻[t in T, p in P] >= 0 # Gas flow entering the pipe
        Q_nm_out⁻[t in T, p in P] >= 0 # Gas flow exiting the pipe
        y_nm[t in T, p in P], Bin
        Qd_cur_gas[t in T, d in D] # load shedding at node n at timestep t [kg]
        LP[t in T, p in P] # Linepack mass
    end
    
    if isempty(compressors)==false # if there are compressors in the network
        @variables m begin    
            Q_c[t in T, c in C] # Gas flow in each compressor [kg/s]
        end
    end

    ################### Objective function #################### 

    if ES.obj_type == "linear"
        @objective(m, Min, dt/3600 *(
            sum(Q_w[t,s]*sources[s].c1 for s in S, t in T) +
            sum(Qd_cur_gas[t,d]*C_cur_load_gas for d in D, t in T) +
            (isempty(compressors) ? 0.0 : sum(compressors[c].compr_cost*(
                pr[t,compressors[c].stop]-pr[t,compressors[c].start]) for c in C, t in T))
            ))
    elseif ES.obj_type == "quadratic"
        @objective(m, Min, dt/3600 *(
            sum(Q_w[t,s]*sources[s].c1 + Q_w[t,s]^2*sources[s].c2 for s in S, t in T) +
            sum(Qd_cur_gas[t,d]*C_cur_load_gas for d in D, t in T) +
            (isempty(compressors) ? 0.0 : sum(compressors[c].compr_cost*(
                pr[t,compressors[c].stop]-pr[t,compressors[c].start]) for c in C, t in T))
            ))
    end

    ################### Constraints ####################
    @constraints m begin
        supply_limits[t in T, s in S],
            sources[s].Q_min <= Q_w[t,s] <= sources[s].Q_max
        pressure_limits[t in T, n in N],
            nodes[n].pr_min <= pr[t,n] <= nodes[n].pr_max 
        bidirectional_flows[t in T, p in P],
            Q_nm[t,p] == Q_nm⁺[t,p] - Q_nm⁻[t,p]
        flow_bounds_pos_upper[t in T, p in P],
            Q_nm⁺[t,p] <= M * y_nm[t,p]
        flow_bounds_neg_upper[t in T, p in P],
            Q_nm⁻[t,p] <= M * (1-y_nm[t,p])
        avg_flow_def_pos[t in T, p in P],
            Q_nm⁺[t,p] == (Q_nm_in⁺[t,p] + Q_nm_out⁺[t,p])/2
        avg_flow_def_neg[t in T, p in P],
            Q_nm⁻[t,p] == (Q_nm_in⁻[t,p] + Q_nm_out⁻[t,p])/2
        linepack_def[t in T, p in P],
            LP[t,p] == pipes[p].S*ϕ⁺[t,p]/2
        linepack_intertemporal_initial[t = [1], p in P],
            Q_nm_in⁺[t,p] - Q_nm_out⁺[t,p] + Q_nm_in⁻[t,p] - Q_nm_out⁻[t,p] == 0
        linepack_intertemporal[t in T, p in P; t > 1],
            LP[t,p] == LP[t-1,p] +
                dt * (Q_nm_in⁺[t,p] - Q_nm_out⁺[t,p] + Q_nm_in⁻[t,p] - Q_nm_out⁻[t,p])
        linepack_end_condition[t = [T[end]], p in P],
            LP[t,p] >= LP[1,p] # pipes[p].LP0
        flow_limits_lower[t in T, p in P],
            -M * (1-y_nm[t,p]) <= Q_nm[t,p]
        flow_limits_upper[t in T, p in P],
            Q_nm[t,p] <= M * y_nm[t,p]
        ϕ⁺_def[t in T, p in P],
            ϕ⁺[t,p]==pr[t,pipes[p].start]+pr[t,pipes[p].stop]
        ϕ⁻_de[t in T, p in P],
            ϕ⁻[t,p]==pr[t,pipes[p].start]-pr[t,pipes[p].stop]

        weymouth_soc_pos[t in T, p in P],
            [pipes[p].Knm^2 * ϕ⁺[t,p] * ϕ⁻[t,p] + M^2 * (1-y_nm[t,p]), Q_nm[t,p]] in SecondOrderCone()
        weymouth_soc_neg[t in T, p in P],
            [pipes[p].Knm^2 * -ϕ⁺[t,p] * ϕ⁻[t,p] + M^2 * y_nm[t,p], Q_nm[t,p]] in SecondOrderCone()
        curtailment_limit[t in T, d in D],
            0 <= Qd_cur_gas[t,d] <= demands[d].Qd_gas[t]
    end


    @expression(m, expr_node_balance[t in T, n in N],
        sum(Q_w[t,s] for s in S if sources[s].location == n) -
        sum(demands[d].Qd_gas[t] for d in D if demands[d].node_load == n) +
        sum(Qd_cur_gas[t,d] for d in D if demands[d].node_load == n) +
        sum(Q_nm_out⁺[t,p] for p in P if pipes[p].stop == n) -
        sum(Q_nm_in⁺[t,p] for p in P if pipes[p].start == n) +
        sum(Q_nm_out⁻[t,p] for p in P if pipes[p].start == n) -
        sum(Q_nm_in⁻[t,p] for p in P if pipes[p].stop == n))

    # Nodal balance equation with and without compressors
    if isempty(compressors) == true # if there are no compressors in the nework
        @constraint(m, node_balance[t in T, n in N], expr_node_balance[t,n] == 0)
    else 
        @constraints m begin
            node_balance[t in T, n in N],
                expr_node_balance[t,n] +
                sum(Q_c[t,c] for c in C if compressors[c].stop == n) -
                sum(Q_c[t,c] for c in C if compressors[c].start == n ) == 0
            compressor_flow[t in T, c in C],
                0 <= Q_c[t,c]
            compressor_limits_l[t in T, c in C],
                compressors[c].CR_min*pr[t,compressors[c].start] <= pr[t,compressors[c].stop] 
            compressor_limits_u[t in T, c in C],
                pr[t,compressors[c].stop] <= compressors[c].CR_max*pr[t,compressors[c].start]
        end
    end

    if ES.print_model == 1
        print(m)
    end

    return m
end


# ####################### Piecewise linear incremental #########################
function build_PWL_incremental_model(ES::GasSystem, no_V)

    nodes, sources, demands, pipes, compressors,
        N, S, D, P, T, C, dt, C_cur_load_gas = get_gas_system_attributes(ES)

    if no_V < 2
        throw(ValueError("Number of set points must be greater than 2."))
    end

    # Derive PWL parameters
    no_K = no_V-1 # number of segments for linearization
    K = 1:no_K 

    ########## Define flow and pressure setpoints ##########
    Q_set = [collect(LinRange(ES.pipes[p].Q_nm_min, ES.pipes[p].Q_nm_max, no_V)) for t in T, p in P]
    pr_set = [collect(LinRange(ES.nodes[n].pr_min, ES.nodes[n].pr_max, no_V)) for t in T, n in N]

    m = Model()
    m = set_solver(m, ES)
    
    #################### Variables ####################
    @variables m begin
        Q_w[t in T, s in S] # Scheduled gas source flowrate [kg] during each timestep
        pr[t in T, n in N] # Nodal pressure  [MPa]
        pr_sq_appr[t in T, n in N] # approximated squared pressure
        y_pr[t in T, n in N, k in K], Bin # auxilary binary variable for pressure linearization
        0 <= δ_pr[t in T, n in N, k in K] <= 1 # auxilary continuous filling variable for pressure linearization
        Q_nm_sq_appr[t in T, p in P] # approximated squared flow
        y_Q[t in T, p in P, k in K], Bin # auxilary binary variable for flow linearization
        0 <= δ_Q[t in T, p in P, k in K] <= 1 # auxilary continuous filling variable for flow linearization
        Q_nm[t in T, p in P] # bidirectional gas flow
        Q_nm_in[t in T, p in P] # Gas flow entering the pipe
        Q_nm_out[t in T, p in P] # Gas flow exiting the pipe
        Qd_cur_gas[t in T, d in D] # load shedding at node n at timestep t [kg]
        LP[t in T, p in P] # Linepack mass
    end
    
    if isempty(compressors)==false # if there are compressors in the nework
        @variables m begin    
        Q_c[t in T, c in C] # Gas flow in each compressor [kg/s]
        end
    end

    ################### Objective function #################### 
    if ES.obj_type == "linear"
        @objective(m, Min, dt/3600 *(
            sum(Q_w[t,s]*sources[s].c1 for s in S, t in T) +
            sum(Qd_cur_gas[t,d]*C_cur_load_gas for d in D, t in T) +
            (isempty(compressors) ? 0.0 : sum(compressors[c].compr_cost*(
                pr[t,compressors[c].stop]-pr[t,compressors[c].start]) for c in C, t in T))
            ))
    elseif ES.obj_type == "quadratic"
        @objective(m, Min, dt/3600 *(
            sum(Q_w[t,s]*sources[s].c1 + Q_w[t,s]^2*sources[s].c2 for s in S, t in T) +
            sum(Qd_cur_gas[t,d]*C_cur_load_gas for d in D, t in T) +
            (isempty(compressors) ? 0.0 : sum(compressors[c].compr_cost*(
                pr[t,compressors[c].stop]-pr[t,compressors[c].start]) for c in C, t in T))
            ))
    end

    ################### Constraints ####################

    @constraints m begin
        # supply_limits[t in T, s in S],
        #     sources[s].Q_min <= Q_w[t,s] <= sources[s].Q_max
        supply_limits_high[t in T, s in S],
            Q_w[t,s] <= sources[s].Q_max
        supply_limits_low[t in T, s in S],
            Q_w[t,s] >= sources[s].Q_min
        # pressure_limits[t in T, n in N],
        #     nodes[n].pr_min <= pr[t,n] <= nodes[n].pr_max
        pressure_limits_high[t in T, n in N],
            pr[t,n] <= nodes[n].pr_max
        pressure_limits_low[t in T, n in N],
            pr[t,n] >= nodes[n].pr_min
        avg_flow_def[t in T, p in P],
            Q_nm[t,p] == (Q_nm_in[t,p] + Q_nm_out[t,p])/2
        linepack_def[t in T, p in P],
            LP[t,p] == pipes[p].S*(pr[t,pipes[p].start]+pr[t,pipes[p].stop])/2
        linepack_intertemporal_initial[t = [1], p in P],
            #ϕ⁺[t,p] /2 / LP[t,p] * (
            Q_nm_in[t,p] - Q_nm_out[t,p] == 0 #why no H0 parameter here?
        linepack_intertemporal[t in T, p in P; t > 1],
            LP[t,p] == LP[t-1,p] +
                dt * (Q_nm_in[t,p] - Q_nm_out[t,p])
        linepack_end_condition[t = [T[end]], p in P],
            LP[t,p] >= LP[1,p] # pipes[p].LP0
        linearized_weymouth[t in T, p in P],
            Q_nm_sq_appr[t,p] == pipes[p].Knm^2 * (pr_sq_appr[t,pipes[p].start] - pr_sq_appr[t,pipes[p].stop])
        # curtailment_limit[t in T, d in D],
        #     0 <= Qd_cur_gas[t,d] <= demands[d].Qd_gas[t]
        curtailment_limit_high[t in T, d in D],
            Qd_cur_gas[t,d] <= demands[d].Qd_gas[t]
        curtailment_limit_low[t in T, d in D],
            Qd_cur_gas[t,d] >= 0
    end
    
    @expression(m, expr_node_balance[t in T, n in N],
        sum(Q_w[t,s] for s in S if sources[s].location == n) -
        sum(demands[d].Qd_gas[t] for d in D if demands[d].node_load == n) +
        sum(Qd_cur_gas[t,d] for d in D if demands[d].node_load == n) +
        sum(Q_nm_out[t,p] for p in P if pipes[p].stop == n) -
        sum(Q_nm_in[t,p] for p in P if pipes[p].start == n))

    # Nodal balance equation with and without compressors
    if isempty(compressors) == true # if there are no compressors in the nework
        @constraint(m, node_balance[t in T, n in N], expr_node_balance[t,n] == 0)
    else 
        @constraints m begin
            node_balance[t in T, n in N],
                expr_node_balance[t,n] +
                sum(Q_c[t,c] for c in C if compressors[c].stop == n) -
                sum(Q_c[t,c] for c in C if compressors[c].start == n ) == 0
            compressor_flow[t in T, c in C],
                0 <= Q_c[t,c]
            compressor_limits_l[t in T, c in C],
                compressors[c].CR_min*pr[t,compressors[c].start] <= pr[t,compressors[c].stop] 
            compressor_limits_u[t in T, c in C],
                pr[t,compressors[c].stop] <= compressors[c].CR_max*pr[t,compressors[c].start]
        end
    end

    ################### Constraints linearization ####################
    @constraints m begin
        Q_nm_sq_appr_def[t in T, p in P],
            Q_nm_sq_appr[t,p] ==  Q_set[t,p][1]*abs(Q_set[t,p][1]) +
                sum((Q_set[t,p][k+1]*abs(Q_set[t,p][k+1])-Q_set[t,p][k]*abs(Q_set[t,p][k]))*δ_Q[t,p,k] for k in K) # This is (15b) , but allows for bidirectionality (abs instead of q^2)
        Q_nm_appr_def[t in T, p in P],
            Q_nm[t,p] == Q_set[t,p][1] + sum((Q_set[t,p][k+1]-Q_set[t,p][k])*δ_Q[t,p,k] for k in K) #This is (15b) in Adrianos formulation 
        δ_Q_limit[t in T, p in P, k in K[1:end-1]],
            δ_Q[t,p,k+1] <= y_Q[t,p,k]
        y_Q_limit[t in T, p in P, k in K[1:end-1]],
            y_Q[t,p,k] <= δ_Q[t,p,k] 
        pr_sq_appr_def[t in T, n in N],
            pr_sq_appr[t,n] ==  pr_set[t,n][1]*abs(pr_set[t,n][1]) +
                sum((pr_set[t,n][k+1]*abs(pr_set[t,n][k+1])-pr_set[t,n][k]*abs(pr_set[t,n][k]))*δ_pr[t,n,k] for k in K) #this is not in adrianos formulation... 
        pr_appr_def[t in T, n in N],
            pr[t,n] == pr_set[t,n][1] + sum((pr_set[t,n][k+1]-pr_set[t,n][k])*δ_pr[t,n,k] for k in K) #this is not in adrianos formulation...
        δ_pr_limit[t in T, n in N, k in K[1:end-1]],
            δ_pr[t,n,k+1] <= y_pr[t,n,k]
        y_pr_limit[t in T, n in N, k in K[1:end-1]],
            y_pr[t,n,k] <= δ_pr[t,n,k] 
        end

    if ES.print_model == 1
        print(m)
    end

    return m

end


####################### MISOCP McCormick #########################
function build_MISOCP_McCormick_model(ES)

    nodes, sources, demands, pipes, compressors,
        N, S, D, P, T, C, dt, C_cur_load_gas = get_gas_system_attributes(ES)
    M = ES.M

    m = Model()
    m = set_solver(m, ES)
 
    #################### Variables ####################
    @variables m begin
        Q_w[t in T, s in S] # Scheduled gas source flowrate [kg] during each timestep
        pr[t in T, n in N] # Nodal pressure  [MPa]
        ϕ⁺[t in T, p in P] # auxiliary variables for pressures ϕ⁺ = pr_m + pr_n
        ϕ⁻[t in T, p in P] # auxiliary variables for pressures ϕ⁻ = pr_m - pr_n
        ψ[t in T, p in P] # auxilary variable used for McCormick envelops, ψ = ϕ⁺ * ϕ⁻
        Q_nm[t in T, p in P] # bidirectional gas flow
        Q_nm⁺[t in T, p in P] >= 0 # Gas flow in each pipe [kg] during each timestep
        Q_nm_in⁺[t in T, p in P] >= 0 # Gas flow entering the pipe
        Q_nm_out⁺[t in T, p in P] >= 0 # Gas flow exiting the pipe
        Q_nm⁻[t in T, p in P] >= 0 # Gas flow in each pipe [kg] during each timestep
        Q_nm_in⁻[t in T, p in P] >= 0 # Gas flow entering the pipe
        Q_nm_out⁻[t in T, p in P] >= 0 # Gas flow exiting the pipe
        y_nm[t in T, p in P], Bin
        Qd_cur_gas[t in T, d in D] # load shedding at node n at timestep t [kg]
        LP[t in T, p in P] # Linepack mass
    end

    
    if isempty(compressors)==false # if there are compressors in the nework
        @variables m begin    
            Q_c[t in T, c in C] # Gas flow in each compressor [kg/s]
        end
    end

    ################### Objective function #################### 
    if ES.obj_type == "linear"
        @objective(m, Min, dt/3600 *(
            sum(Q_w[t,s]*sources[s].c1 for s in S, t in T) +
            sum(Qd_cur_gas[t,d]*C_cur_load_gas for d in D, t in T) +
            (isempty(compressors) ? 0.0 : sum(compressors[c].compr_cost*(
                pr[t,compressors[c].stop]-pr[t,compressors[c].start]) for c in C, t in T))
            ))
    elseif ES.obj_type == "quadratic"
        @objective(m, Min, dt/3600 *(
            sum(Q_w[t,s]*sources[s].c1 + Q_w[t,s]^2*sources[s].c2 for s in S, t in T) +
            sum(Qd_cur_gas[t,d]*C_cur_load_gas for d in D, t in T) +
            (isempty(compressors) ? 0.0 : sum(compressors[c].compr_cost*(
                pr[t,compressors[c].stop]-pr[t,compressors[c].start]) for c in C, t in T))
            ))
    end

    ################### Constraints ####################
    @constraints m begin
        supply_limits[t in T, s in S],
            sources[s].Q_min <= Q_w[t,s] <= sources[s].Q_max

        pressure_limits_aux_pos[t in T, p in P],
            pipes[p].ϕ⁺_min <= ϕ⁺[t,p] <= pipes[p].ϕ⁺_max
        pressure_limits_aux_neg[t in T, p in P],
            pipes[p].ϕ⁻_min <= ϕ⁻[t,p] <= pipes[p].ϕ⁻_max


        bidirectional_flows[t in T, p in P],
            Q_nm[t,p] == Q_nm⁺[t,p] - Q_nm⁻[t,p]
        flow_bounds_pos_upper[t in T, p in P],
            Q_nm⁺[t,p] <= M * y_nm[t,p]
        flow_bounds_neg_upper[t in T, p in P],
            Q_nm⁻[t,p] <= M * (1-y_nm[t,p])
        avg_flow_def_pos[t in T, p in P],
            Q_nm⁺[t,p] == (Q_nm_in⁺[t,p] + Q_nm_out⁺[t,p])/2
        avg_flow_def_neg[t in T, p in P],
            Q_nm⁻[t,p] == (Q_nm_in⁻[t,p] + Q_nm_out⁻[t,p])/2

        linepack_def[t in T, p in P],
            LP[t,p] == pipes[p].S*ϕ⁺[t,p]/2

        linepack_intertemporal_initial[t = [1], p in P],
            Q_nm_in⁺[t,p] - Q_nm_out⁺[t,p] + Q_nm_in⁻[t,p] - Q_nm_out⁻[t,p] == 0

        linepack_intertemporal[t in T, p in P; t > 1],
            LP[t,p] == LP[t-1,p] +
                dt * (Q_nm_in⁺[t,p] - Q_nm_out⁺[t,p] + Q_nm_in⁻[t,p] - Q_nm_out⁻[t,p])

        # linepack_end_condition[t = [T[end]], p in P],
        #     LP[t,p] >= LP[1,p] # pipes[p].LP0

        linepack_end_condition, # Restore linepack at the end of the time horizon
            sum(Q_w[t,s] for t in T, s in S) >= sum(demands[d].Qd_gas[t] for t in T, d in D)

        flow_limits_lower[t in T, p in P],
            -M * (1-y_nm[t,p]) <= Q_nm[t,p]
        flow_limits_upper[t in T, p in P],
            Q_nm[t,p] <= M * y_nm[t,p]



        ϕ_limits_1[t in T, p in P],
            nodes[pipes[p].start].pr_min <= 0.5*(ϕ⁺[t,p]+ϕ⁻[t,p]) <= nodes[pipes[p].start].pr_max
        ϕ_limits_2[t in T, p in P],
            nodes[pipes[p].stop].pr_min <= 0.5*(ϕ⁺[t,p]-ϕ⁻[t,p]) <= nodes[pipes[p].stop].pr_max

        ϕ⁺_def[t in T, p in P],
            ϕ⁺[t,p]==pr[t,pipes[p].start]+pr[t,pipes[p].stop]
        ϕ⁻_def[t in T, p in P],
            ϕ⁻[t,p]==pr[t,pipes[p].start]-pr[t,pipes[p].stop]

        weymouth_soc_pos[t in T, p in P],
            [pipes[p].Knm^2 * ψ[t,p] + M^2 * (1-y_nm[t,p]), Q_nm[t,p]] in SecondOrderCone()
        weymouth_soc_neg[t in T, p in P],
            [-pipes[p].Knm^2 * ψ[t,p] + M^2 * y_nm[t,p], Q_nm[t,p]] in SecondOrderCone()
        McCormick_1[t in T, p in P],
            ψ[t,p] >= pipes[p].ϕ⁺_min * ϕ⁻[t,p] + ϕ⁺[t,p] * pipes[p].ϕ⁻_min - pipes[p].ϕ⁺_min * pipes[p].ϕ⁻_min
        McCormick_2[t in T, p in P],
            ψ[t,p] >= pipes[p].ϕ⁺_max * ϕ⁻[t,p] + ϕ⁺[t,p] * pipes[p].ϕ⁻_max - pipes[p].ϕ⁺_max * pipes[p].ϕ⁻_max
        McCormick_3[t in T, p in P],
            ψ[t,p] <= pipes[p].ϕ⁺_max * ϕ⁻[t,p] + ϕ⁺[t,p] * pipes[p].ϕ⁻_min - pipes[p].ϕ⁺_max * pipes[p].ϕ⁻_min
        McCormick_4[t in T, p in P],
            ψ[t,p] >= pipes[p].ϕ⁺_min * ϕ⁻[t,p] + ϕ⁺[t,p] * pipes[p].ϕ⁻_max - pipes[p].ϕ⁺_min * pipes[p].ϕ⁻_max

        curtailment_limit[t in T, d in D],
            0 <= Qd_cur_gas[t,d] <= demands[d].Qd_gas[t]
    end


    @expression(m, expr_node_balance[t in T, n in N],
        sum(Q_w[t,s] for s in S if sources[s].location == n) -
        sum(demands[d].Qd_gas[t] for d in D if demands[d].node_load == n) +
        sum(Qd_cur_gas[t,d] for d in D if demands[d].node_load == n) +
        sum(Q_nm_out⁺[t,p] for p in P if pipes[p].stop == n) -
        sum(Q_nm_in⁺[t,p] for p in P if pipes[p].start == n) +
        sum(Q_nm_out⁻[t,p] for p in P if pipes[p].start == n) -
        sum(Q_nm_in⁻[t,p] for p in P if pipes[p].stop == n))

    # Nodal balance equation with and without compressors
    if isempty(compressors) == true # if there are no compressors in the nework
        @constraint(m, node_balance[t in T, n in N], expr_node_balance[t,n] == 0)
    else 
        @constraints m begin
            node_balance[t in T, n in N],
                expr_node_balance[t,n] +
                sum(Q_c[t,c] for c in C if compressors[c].stop == n) -
                sum(Q_c[t,c] for c in C if compressors[c].start == n ) == 0
            compressor_flow[t in T, c in C],
                0 <= Q_c[t,c]
            compressor_limits_l[t in T, c in C],
                compressors[c].CR_min*pr[t,compressors[c].start] <= pr[t,compressors[c].stop] 
            compressor_limits_u[t in T, c in C],
                pr[t,compressors[c].stop] <= compressors[c].CR_max*pr[t,compressors[c].start]
        end
    end

    if ES.print_model == 1
        print(m)
    end

    return m 
end

function build_NN_constrained_model(ES, NN; M_up = 100, M_low = -100, individual_pipes = false, heuristic_level = 0.05)

    nodes, sources, demands, pipes, compressors,
        N, S, D, P, T, C, dt, C_cur_load_gas = get_gas_system_attributes(ES)

    m = Model()
    # m.setParam('BranchDir', 1)

    # m=direct_model(Gurobi.Optimizer())
    
    m = set_solver(m, ES)
    if solver_name(m) == "Gurobi"
        # set_optimizer_attribute(m, "BranchDir", 1)
        # set_optimizer_attribute(m, "Heuristics", heuristic_level)
        # set_optimizer_attribute(m, "MIPFocus", 3)
    end
    if solver_name(m) == "CPLEX"
        # set_optimizer_attribute(m, "CPXPARAM_Emphasis_MIP", 1) #Emphasise feasible solutions
        # set_optimizer_attribute(m, "Heuristics", heuristic_level)
    end
    # set_optimizer_attribute(m, "MIPFocus", 3)
    # set_optimizer_attribute(m, "Heuristics", heuristic_level)
    #if isnothing(individual_pipes)
    if individual_pipes
        ind_pipe_dir = joinpath(dirname(@__DIR__), "inputs", "neural_network_par","parameters", "ipopt", "individual_pipes_run_2023_05_07_15_05_09")
        ind_pipe_files = readdir(ind_pipe_dir)
        
        W = Any[]
        B = Any[]
        K = Any[]
        α = Any[]
        y_hat_min = Any[]
        y_hat_max = Any[]
        prn_mean = Any[]
        prn_std = Any[]
        prm_mean = Any[]
        prm_std = Any[]
        Q_mean = Any[]
        Q_std = Any[]
        #maybe there is a more clean way to do this, than just doing lists... but it works
        println(ind_pipe_files)
        for (i,x) in enumerate(ind_pipe_files)
            pipe = "pipe$i.json"
            global file
            file = ""
            for p in ind_pipe_files #loop over all files to extract the correct one... 
                if occursin(pipe, p)
                    println(pipe, p)
                    global file
                    file = p
                    # println(file)
                end
            end
            inp_file = joinpath(ind_pipe_dir, file)
            println("\npipe$i ;", inp_file)
            NN_data = JSON.parsefile(inp_file)
            Wp, Bp, Kp, αp, y_hat_min_inp, y_hat_max_inp, mn, sd = read_NN_parameters(NN_data)
            #combine all output for each pipe into a single dictionary
            Ne = [1:size(Wp[2])[2],1:size(Wp[3])[2]] #number of neurons in each layer
            NeI = 1:size(Wp[1])[2] #number of neurons in the input layer

            #bounds on the input to the NN
            y_hat_minp = min.(y_hat_min_inp,-0.0001)#.*1.05#minimum value of y (this should be a parameter of get_gas_system_attributes)
            y_hat_maxp = max.(y_hat_max_inp,0.0001)#.*1.05 #maximum value of y (this should be a parameter of get_gas_system_attributes)
            y_hat_minp = max.(y_hat_minp, M_low)
            y_hat_maxp = min.(y_hat_maxp, M_up)
            # println("y_hat_min $y_hat_min")
            # println("y_hat_maxp $y_hat_maxp")

            #values for scaling
            #pressure at node n (input feature 1)
            prn_meanp = mn[1] #mean value of pr (this should be a parameter of get_gas_system_attributes)
            prn_stdp = sd[1] #standard deviation of pr (this should be a parameter of get_gas_system_attributes)
            #pressure at node m (input feature 2)
            prm_meanp = mn[2] #mean value of pr (this should be a parameter of get_gas_system_attributes)
            prm_stdp = sd[2] #standard deviation of pr (this should be a parameter of get_gas_system_attributes)
            #flowrate Q_nm (target feature)
            Q_meanp = mn[3] #standard deviation of Q (this should be a parameter of get_gas_system_attributes)
            Q_stdp = sd[3] #mean value of Q (this should be a parameter of get_gas_system_attributes)
            # nn_pipe_dict[i] = (W, B, K, α, y_hat_min_in, y_hat_max_in, mn, sd)
            push!(W, Wp);
            push!(B, Bp);
            push!(K, Kp);
            push!(α, αp);
            push!(y_hat_min, y_hat_minp);
            push!(y_hat_max, y_hat_maxp);
            push!(prn_mean, prn_meanp);
            push!(prn_std, prn_stdp);
            push!(prm_mean, prm_meanp);
            push!(prm_std, prm_stdp);
            push!(Q_mean, Q_meanp);
            push!(Q_std, Q_stdp);    
            #should add NeI as well... 
        end
        @variables m begin
            #Neural network variables:
            #define for each timestep, each pipe, each layer, each neuron
            y0[t=1:length(T), p=1:length(P), ne=1:length(NeI)]  #initial input (needs to be scaled from pr) #Maybe add this as middle step for simplicity
            y[t=1:length(T), p=1:length(P), k=1:length(K[1]), ne=1:length(Ne[1])] #input for next neuron
            y_hat[t=1:length(T), p=1:length(P), k=1:length(K[1]), ne=1:length(Ne[1])] #output of each neuron
            z[t=1:length(T), p=1:length(P)] #final output of NN (this needs to be converted to Q_nm)
            xn[t=1:length(T), p=1:length(P), k=1:length(K[1]), ne=1:length(Ne[1])] <= 0 #negative part of input
            xp[t=1:length(T), p=1:length(P), k=1:length(K[1]), ne=1:length(Ne[1])] >= 0 #positive part of input
            b[t=1:length(T), p=1:length(P), k=1:length(K[1]), ne=1:length(Ne[1])], Bin 
        end

        # y_hat[t,p,k+1,:] .== W[p][k+1]*y[t,p,k,:]+B[p][k+1]
    else
        # NN hyperparameters, Weights and biases:
        println("nnpath $NN")
        NN_data = JSON.parsefile(NN) #change this to be an input parameter... 

        # W, B, 1:K-1, α, y1_bounds, y2_bounds, mn, std
        W, B, K, α, y_hat_min_in, y_hat_max_in, mn, sd = read_NN_parameters(NN_data) #read the weights and biases from a dictionary
        # W = (ones(5,2),ones(5,5),ones(1,5)) #matrix of shape([[y[k] Ne[k]], K+1]) weights of each neuron in each layer (this should be a parameter of get_gas_system_attributes)
        # B = (ones(5),ones(5),ones(1)) #matrix of biases shape([[Ne[k]], K+1]) bias of each neuron in each layer (this should be a parameter of get_gas_system_attributes)
        Ne = [1:size(W[2])[2],1:size(W[2])[2]] #number of neurons in each layer
        NeI = 1:size(W[1])[2] #number of neurons in the input layer

        #bounds on the input to the NN
        y_hat_min = min.(y_hat_min_in,-0.0001)#.*1.05#minimum value of y (this should be a parameter of get_gas_system_attributes)
        y_hat_max = max.(y_hat_max_in,0.0001)#.*1.05 #maximum value of y (this should be a parameter of get_gas_system_attributes)
        y_hat_min = max.(y_hat_min, M_low)
        y_hat_max = min.(y_hat_max, M_up)
        println("y_hat_min $y_hat_min")
        println("y_hat_max $y_hat_max")

        #values for scaling
        #pressure at node n (input feature 1)
        prn_mean = mn[1] #mean value of pr (this should be a parameter of get_gas_system_attributes)
        prn_std = sd[1] #standard deviation of pr (this should be a parameter of get_gas_system_attributes)
        #pressure at node m (input feature 2)
        prm_mean = mn[2] #mean value of pr (this should be a parameter of get_gas_system_attributes)
        prm_std = sd[2] #standard deviation of pr (this should be a parameter of get_gas_system_attributes)
        #flowrate Q_nm (target feature)
        Q_mean = mn[3] #standard deviation of Q (this should be a parameter of get_gas_system_attributes)
        Q_std = sd[3] #mean value of Q (this should be a parameter of get_gas_system_attributes)

        @variables m begin
            #Neural network variables:
            #define for each timestep, each pipe, each layer, each neuron
            y0[t=1:length(T), p=1:length(P), ne=1:length(NeI)]  #initial input (needs to be scaled from pr) #Maybe add this as middle step for simplicity
            y[t=1:length(T), p=1:length(P), k=1:length(K), ne=1:length(Ne[1])] #input for next neuron
            y_hat[t=1:length(T), p=1:length(P), k=1:length(K), ne=1:length(Ne[1])] #output of each neuron
            z[t=1:length(T), p=1:length(P)] #final output of NN (this needs to be converted to Q_nm)
            xn[t=1:length(T), p=1:length(P), k=1:length(K), ne=1:length(Ne[1])] <= 0 #negative part of input
            xp[t=1:length(T), p=1:length(P), k=1:length(K), ne=1:length(Ne[1])] >= 0 #positive part of input
            b[t=1:length(T), p=1:length(P), k=1:length(K), ne=1:length(Ne[1])], Bin 
        end

    end

    #################### Variables ####################
    @variables m begin #it makes a difference to say t in T or just T...
        Q_w[t in T, s in S] # Scheduled gas source flowrate [kg] during each timestep
        pr[t in T, n in N] # Nodal pressure  [MPa] #this also means that compressors are a subset of pipes... 
        # pr_sq_appr[t in T, n in N] # approximated squared pressure
        # y_pr[t in T, n in N, k in K], Bin # auxilary binary variable for pressure linearization
        # 0 <= δ_pr[t in T, n in N, k in K] <= 1 # auxilary continuous filling variable for pressure linearization
        # Q_nm_sq_appr[t in T, p in P] # approximated squared flow
        # y_Q[t in T, p in P, k in K], Bin # auxilary binary variable for flow linearization
        # 0 <= δ_Q[t in T, p in P, k in K] <= 1 # auxilary continuous filling variable for flow linearization
        Q_nm[t in T, p in P] # bidirectional gas flow
        Q_nm_in[t in T,p in P] # Gas flow entering the pipe
        Q_nm_out[t in T, p in P] # Gas flow exiting the pipe
        Qd_cur_gas[t in T, d in D] # load shedding at node n at timestep t [kg]
        LP[t in T, p in P] # Linepack mass
    end

    if isempty(compressors)==false # if there are compressors in the nework
        @variable(m, Q_c[t in T, c in C] ) # Gas flow in each compressor [kg/s]
    end

    #################### Objective function ####################

    if ES.obj_type == "linear"
        @objective(m, Min, dt/3600 *(
            sum(Q_w[t,s]*sources[s].c1 for s in S, t in T) +
            sum(Qd_cur_gas[t,d]*C_cur_load_gas for d in D, t in T) +
            (isempty(compressors) ? 0.0 : sum(compressors[c].compr_cost*(
                pr[t,compressors[c].stop]-pr[t,compressors[c].start]) for c in C, t in T))
            ))
    elseif ES.obj_type == "quadratic"
        @objective(m, Min, dt/3600 *(
            sum(Q_w[t,s]*sources[s].c1 + Q_w[t,s]^2*sources[s].c2 for s in S, t in T) +
            sum(Qd_cur_gas[t,d]*C_cur_load_gas for d in D, t in T) +
            (isempty(compressors) ? 0.0 : sum(compressors[c].compr_cost*(
                pr[t,compressors[c].stop]-pr[t,compressors[c].start]) for c in C, t in T))
            ))
    end


    #################### Constraints ####################

    @constraints m begin
        # supply_limits[t in T, s in S],
        #     sources[s].Q_min <= Q_w[t,s] <= sources[s].Q_max
        supply_limits_low[t in T, s in S],
            sources[s].Q_min <= Q_w[t,s]
        supply_limits_high[t in T, s in S],
            Q_w[t,s] <= sources[s].Q_max
        # pressure_limits[t in T, n in N],
        #     nodes[n].pr_min <= pr[t,n] <= nodes[n].pr_max
        pressure_limits_low[t in T, n in N],
            nodes[n].pr_min <= pr[t,n]
        pressure_limits_high[t in T, n in N],
            pr[t,n] <= nodes[n].pr_max
        avg_flow_def[t in T, p in P],
            Q_nm[t,p] == (Q_nm_in[t,p] + Q_nm_out[t,p])/2
        linepack_def[t in T, p in P],
            LP[t,p] == pipes[p].S*(pr[t,pipes[p].start]+pr[t,pipes[p].stop])/2
        linepack_intertemporal_initial[t = [1], p in P], #why the brackets around 1?
            #ϕ⁺[t,p] /2 / LP[t,p] * (
            Q_nm_in[t,p] - Q_nm_out[t,p] == 0 #Why no H0 parameter here?
        linepack_intertemporal[t in T, p in P; t > 1],
            LP[t,p] == LP[t-1,p] +
                dt * (Q_nm_in[t,p] - Q_nm_out[t,p])
        linepack_end_condition[t = [T[end]], p in P],
            LP[t,p] >= LP[1,p] # pipes[p].LP0
        # linearized_weymouth[t in T, p in P], #replace with Neural Network constraints
        #     Q_nm_sq_appr[t,p] == pipes[p].Knm^2 * (pr_sq_appr[t,pipes[p].start] - pr_sq_appr[t,pipes[p].stop]) #replace with Neural Network constraints
        # curtailment_limit[t in T, d in D],
        #     0 <= Qd_cur_gas[t,d] <= demands[d].Qd_gas[t]
        curtailment_limit_low[t in T, d in D],
            0 <= Qd_cur_gas[t,d]
        curtailment_limit_high[t in T, d in D],
            Qd_cur_gas[t,d] <= demands[d].Qd_gas[t]
    end

    @expression(m, expr_node_balance[t in T, n in N],
        sum(Q_w[t,s] for s in S if sources[s].location == n) -
        sum(demands[d].Qd_gas[t] for d in D if demands[d].node_load == n) +
        sum(Qd_cur_gas[t,d] for d in D if demands[d].node_load == n) +
        sum(Q_nm_out[t,p] for p in P if pipes[p].stop == n) -
        sum(Q_nm_in[t,p] for p in P if pipes[p].start == n))

    # Nodal balance equation with and without compressors
    if isempty(compressors) == true # if there are no compressors in the nework
        @constraint(m, node_balance[t in T, n in N], expr_node_balance[t,n] == 0)
    else 
        @constraints m begin
            node_balance[t in T, n in N],
                expr_node_balance[t,n] +
                sum(Q_c[t,c] for c in C if compressors[c].stop == n) -
                sum(Q_c[t,c] for c in C if compressors[c].start == n ) == 0 #Q^C_[t,n,m]
            compressor_flow[t in T, c in C],
                0 <= Q_c[t,c]
            compressor_limits_l[t in T, c in C],
                compressors[c].CR_min*pr[t,compressors[c].start] <= pr[t,compressors[c].stop] 
            compressor_limits_u[t in T, c in C],
                pr[t,compressors[c].stop] <= compressors[c].CR_max*pr[t,compressors[c].start]
        end
    end


        ################### Constraints Neural Network linearization ####################
    if individual_pipes
        """
        Using a neural network for each pipe
        
            NOTE: if the neural networks does not have the same dimensions we will have a problem... this is because the variable needs to be matrix{variableRef} type. The code structure is not set up to handle different neural networks structures.
            If different number of layers or neurons is used for each pipe, then we need to use a different approach. We might need to hardcode the variables for each pipe in that case... or somethings. 
        """
        println("using a neural network for each pipe")
        @constraints m begin
            #This is the Neural Network approximation of the Weymouth equation
            #NOTE: dot product (inner product) is computed using * (make sure that dimensions align...)
            #input layer
            input_var[t in T, p in P],
                y0[t,p,:] .== [(pr[t,pipes[p].start]-prn_mean[p])/prn_std[p],(pr[t,pipes[p].stop]-prm_mean[p])/prm_std[p]] #check mean/std calcs. might need to transpose...
            input_layer[t in T, p in P],
                y_hat[t,p,1,:] .== W[p][1]*y0[t,p,:]+B[p][1]
            #hidden layer
            hidden_layers[t in T, p in P, k in 1:length(K[1])-1], #better way to write 1:K-1?
                y_hat[t,p,k+1,:] .== W[p][k+1]*y[t,p,k,:]+B[p][k+1]

            #activation function constraints
            neg_lim[t in T, p in P, k in K[1]], #check if i have to define N in constraint def... 
                y_hat_min[p][:,k] .*b[t,p,k,:] .<= xn[t,p,k,:] 
            pos_lim[t in T, p in P, k in K[1]],
                xp[t,p,k,:] .<= y_hat_max[p][:,k].*(1 .- b[t,p,k,:]) #y_hat_max...
            yhat_con[t in T, p in P, k in K[1]],
                y_hat[t,p,k,:] .== xn[t,p,k,:] .+ xp[t,p,k,:]
            relu_con[t in T, p in P, k in K[1]],
                y[t,p,k,:] .== α[p].*xn[t,p,k,:].+xp[t,p,k,:]

            #output layer
            output_layer[t in T, p in P, k=length(K[1])],
                z[t,p] .== W[p][k+1]*y[t,p,k,:]+B[p][k+1] #why do i need the .== here?

            #rescale to Q_nm values
            rescale[t in T, p in P],
                Q_nm[t,p] == (z[t,p]*Q_std[p] + Q_mean[p]) # *pipes[p].Knm #dont need to rescale by Knm since we are using the NN approximation of the Weymouth equation for each pipe with k included

            end
    else
        
        @constraints m begin
            #This is the Neural Network approximation of the Weymouth equation
            #NOTE: dot product (inner product) is computed using * (make sure that dimensions align...)
            #input layer
            input_var[t in T, p in P],
                y0[t,p,:] .== [(pr[t,pipes[p].start]-prn_mean)/prn_std,(pr[t,pipes[p].stop]-prm_mean)/prm_std] #check mean/std calcs. might need to transpose...
            input_layer[t in T, p in P],
                y_hat[t,p,1,:] .== W[1]*y0[t,p,:]+B[1]
            #hidden layer
            hidden_layers[t in T, p in P, k in 1:length(K)-1], #better way to write 1:K-1?
                y_hat[t,p,k+1,:] .== W[k+1]*y[t,p,k,:]+B[k+1]

            #activation function constraints
            neg_lim[t in T, p in P, k in K], #check if i have to define N in constraint def... 
                y_hat_min[:,k] .*b[t,p,k,:] .<= xn[t,p,k,:] 
            pos_lim[t in T, p in P, k in K],
                xp[t,p,k,:] .<= y_hat_max[:,k].*(1 .- b[t,p,k,:]) #y_hat_max...
            yhat_con[t in T, p in P, k in K],
                y_hat[t,p,k,:] .== xn[t,p,k,:] .+ xp[t,p,k,:]
            relu_con[t in T, p in P, k in K],
                y[t,p,k,:] .== α.*xn[t,p,k,:].+xp[t,p,k,:]

            #output layer
            output_layer[t in T, p in P, k=length(K)],
                z[t,p] .== W[k+1]*y[t,p,k,:]+B[k+1] #why do i need the .== here?

            #rescale to Q_nm values
            rescale[t in T, p in P],
                Q_nm[t,p] == (z[t,p]*Q_std + Q_mean)*pipes[p].Knm 

            end
    end
    
        if ES.print_model == 1
            print(m)
        end
    
        return m
    
    end