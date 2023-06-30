include("energy_system.jl")

################### PDE MODEL ####################
# Power: DC power flow
# Gas: PDE model

################### PDE MODEL ####################
function build_PDE_model(ES::IntegratedEnergySystem)

    ########## Gas system ##########
    nodes, sources, demands, pipes, compressors,
        N, S, D, P, T, C, dt, C_cur_load_gas = get_gas_system_attributes(ES)

    # PDE Type: static, quasi-dynamic, dynamic (transient)
    ut1, ut2 = PDE_type_parameters(ES.PDE_type)

    ########## Power system ##########
    generators, windgenerators, demands_EL, lines,
        B, G, J, G, D_el, L = get_power_system_attributes(ES)
    GFPPs = [g for g in G if generators[g].type == "NGFPP"]
    nonGFPPs = [g for g in G if generators[g].type == "non-NGFPP"]

    #################### Initializing model ####################
    m = Model()
    set_solver(m, ES)

    #################### Variables ####################
    @variables m begin
        Q_w[t in T, s in S] # Scheduled gas source flowrate [kg/s] during each timestep
        pr[t in T, n in N] # Nodal pressure  [MPa]
        Q_nm[t in T, p in P] # Gas flow in each pipe [kg/s] during each timestep
        Q_nm_in[t in T, p in P] # Gas flow entering the pipe
        Q_nm_out[t in T, p in P] # Gas flow exiting the pipe
        pr_avg[t in T, p in P] # Average pressure in each pipe (MPa)
        Qd_cur_gas[t in T, d in D] # load shedding at node n at timestep t [kg/s]
        Qd_gfpp[t in T, gg in GFPPs] # Gas consumption from GFPPs [kg/s]
        p[g in G, t in T] >= 0 # power dispatched by generators[MW]
        w[j in J, t in T] >= 0 # power dispatched by wind turbines [MW]
    end

    if isempty(compressors)==false # if there are compressors in the nework
        @variable(m, Q_c[t in T, c in C]) # Gas flow in each compressor [kg/s]
    end

    ################### Objective function ####################

    # Cost of gas supply (including gas used for GFPPs) and curtailment
    cost_gas = dt/3600 *(sum(Q_w[t,s]*sources[s].c1 + Q_w[t,s]^2*sources[s].c2 for s in S, t in T) +
        sum(Qd_cur_gas[t,d]*C_cur_load_gas for d in D, t in T)) #+
        # (isempty(compressors) ? 0.0 : sum(compressors[c].compr_cost*(
        #     pr[t,compressors[c].stop]-pr[t,compressors[c].start]) for c in C, t in T))
    # Cost of running non-GFPPs
    cost_el = dt/3600 *(sum(generators[g].c_el*p[g, t] for g in nonGFPPs, t in T))

    @objective(m, Min, cost_gas + cost_el) 

    ################### Constraints ####################
    @constraints m begin
        supply_limits[t in T, s in S],
            sources[s].Q_min <= Q_w[t,s] <= sources[s].Q_max
        # Only the pressure for the original nodes must be constrained
        pressure_limits[t in T, n in N],
            nodes[n].pr_min <= pr[t,n] <= nodes[n].pr_max
        curtailment_limit[t in T, d in D],
            0 <= Qd_cur_gas[t,d] <= demands[d].Qd_gas[t]
    end

    @expression(m, expr_node_balance[t in T, n in N],
        sum(Q_w[t,s] for s in S if sources[s].location == n) -
        sum(demands[d].Qd_gas[t] for d in D if demands[d].node_load == n) +
        sum(Qd_cur_gas[t,d] for d in D if demands[d].node_load == n) +
        sum(Q_nm_out[t,p] for p in P if pipes[p].stop == n) -
        sum(Q_nm_in[t,p] for p in P if pipes[p].start == n) -
        sum(generators[gg].eta*p[gg, t] for gg in GFPPs if generators[gg].NGnode == n))

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

    @constraints m begin
        conservation_momentum_initial[t=[1], p in P],
            pr_avg[t,p]*pipes[p].V_p*(pr[t,pipes[p].stop]-pr[t,pipes[p].start])/pipes[p].LL +
            pipes[p].V_f*(Q_nm[t,p]*Q_nm[t,p]) == 0
        conservation_momentum[t in T, p in P; t >= 2],
            pr_avg[t,p]*ut2*(Q_nm[t,p]-Q_nm[t-1,p])/dt +
            pr_avg[t,p]*pipes[p].V_p*(pr[t,pipes[p].stop]-pr[t,pipes[p].start])/pipes[p].LL +
            pipes[p].V_f*(Q_nm[t,p]*Q_nm[t,p]) == 0
    end

    # Restore linepack at the end of the time horizon
    @constraint(m, linepack_end_condition,
        sum(Q_w[t,s] for t in T, s in S) >= 
            sum(demands[d].Qd_gas[t] for t in T, d in D) +
            sum(generators[gg].eta*p[gg, t] for t in T, gg in GFPPs))
  
    @constraints m begin
        GenerationLimitUp[g in G, t in T],
            p[g, t] <= generators[g].pmax
        # RampUp[g in G, t in T[2:end]], 
        #     p[g,t]-p[g, t-1] <= generators[g].pup/3600*dt
        # RampDown[g in G, t in T[2:end]],
        #     p[g,t]-p[g, t-1] >= -generators[g].pdown/3600*dt
        WindLimitUp[j in J, t in T],
            w[j, t] <= windgenerators[j].P_w[t]

        # Formulation with PTDF matrix
        PowerBalance[t in T],
            sum(p[g, t] for g in G)  + sum(w[j,t] for j in J) - sum(demands_EL[d].D_el[t] for d in D_el) == 0 
        PowerFlowUp[l in L, t in T],
            lines[l].ptdf[[generators[g].ELnode for g in G]]'*([p[g, t] for g in G]) +
            lines[l].ptdf[[windgenerators[j].ELnode for j in J]]'*([w[j, t] for j in J]) -
            lines[l].ptdf[[demands_EL[d].ELnode for d in D_el]]'*([demands_EL[d].D_el[t] for d in D_el]) <= lines[l].fmax
        
        PowerFlowDown[l in L, t in T],
            lines[l].ptdf[[generators[g].ELnode for g in G]]'*([p[g, t] for g in G]) +
            lines[l].ptdf[[windgenerators[j].ELnode for j in J]]'*([w[j, t] for j in J]) -
            lines[l].ptdf[[demands_EL[d].ELnode for d in D_el]]'*([demands_EL[d].D_el[t] for d in D_el]) >= -lines[l].fmax 
    end
    
    if ES.print_model == 1
        print(m)
    end

    return m
end


################### PDE MODEL ####################
function build_PDE_model_SLP(ES::IntegratedEnergySystem)

    ########## Gas system ##########
    nodes, sources, demands, pipes, compressors,
        N, S, D, P, T, C, dt, C_cur_load_gas = get_gas_system_attributes(ES)

    # PDE Type: static, quasi-dynamic, dynamic (transient)
    ut1, ut2 = PDE_type_parameters(ES.PDE_type)

    ########## Power system ##########
    generators, windgenerators, demands_EL, lines,
        B, G, J, G, D_el, L = get_power_system_attributes(ES)
    GFPPs = [g for g in G if generators[g].type == "NGFPP"]
    nonGFPPs = [g for g in G if generators[g].type == "non-NGFPP"]

    #################### Initializing model ####################
    m = Model()
    set_solver(m, ES)

    #################### Variables ####################
    @variables m begin
        Q_w[t in T, s in S] # Scheduled gas source flowrate [kg/s] during each timestep
        pr[t in T, n in N] # Nodal pressure  [MPa]
        Q_nm[t in T, p in P] # Gas flow in each pipe [kg/s] during each timestep
        Q_nm_in[t in T, p in P] # Gas flow entering the pipe
        Q_nm_out[t in T, p in P] # Gas flow exiting the pipe
        pr_avg[t in T, p in P] # Average pressure in each pipe (MPa)
        Qd_cur_gas[t in T, d in D] # load shedding at node n at timestep t [kg/s]
        γ[t in T, p in P] # Auxiliary variable Taylor Series Expansion
        Qd_gfpp[t in T, gg in GFPPs] # Gas consumption from GFPPs [kg/s]
        p[g in G, t in T] >= 0 # power dispatched by generators[MW]
        w[j in J, t in T] >= 0 # power dispatched by wind turbines [MW]
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
        # Only the pressure for the original nodes must be constrained
        pressure_limits[t in T, n in N],
            nodes[n].pr_min <= pr[t,n] <= nodes[n].pr_max
        curtailment_limit[t in T, d in D],
            0 <= Qd_cur_gas[t,d] <= demands[d].Qd_gas[t]
    end

    @expression(m, expr_node_balance[t in T, n in N],
        sum(Q_w[t,s] for s in S if sources[s].location == n) -
        sum(demands[d].Qd_gas[t] for d in D if demands[d].node_load == n) +
        sum(Qd_cur_gas[t,d] for d in D if demands[d].node_load == n) +
        sum(Q_nm_out[t,p] for p in P if pipes[p].stop == n) -
        sum(Q_nm_in[t,p] for p in P if pipes[p].start == n) -
        sum(generators[gg].eta*p[gg, t] for gg in GFPPs if generators[gg].NGnode == n))

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

    @constraints m begin
        conservation_momentum_initial[t=[1], p in P],
            pipes[p].V_p * (pr[t,pipes[p].stop]-pr[t,pipes[p].start])/pipes[p].LL +
            pipes[p].V_f * γ[t,p] == 0
        conservation_momentum[t in T, p in P; t >= 2],
            ut2 * (Q_nm[t,p]-Q_nm[t-1,p])/dt +
            pipes[p].V_p * (pr[t,pipes[p].stop]-pr[t,pipes[p].start])/pipes[p].LL +
            pipes[p].V_f * γ[t,p] == 0
    end

    # Restore linepack at the end of the time horizon
    @constraint(m, linepack_end_condition,
        sum(Q_w[t,s] for t in T, s in S) >= 
            sum(demands[d].Qd_gas[t] for t in T, d in D) +
            sum(generators[gg].eta*p[gg, t] for t in T, gg in GFPPs))
  
    @constraints m begin
        GenerationLimitUp[g in G, t in T],
            p[g, t] <= generators[g].pmax
        # RampUp[g in G, t in T[2:end]], 
        #     p[g,t]-p[g, t-1] <= generators[g].pup/3600*dt
        # RampDown[g in G, t in T[2:end]],
        #     p[g,t]-p[g, t-1] >= -generators[g].pdown/3600*dt
        WindLimitUp[j in J, t in T],
            w[j, t] <= windgenerators[j].P_w[t]

        # Formulation with PTDF matrix
        PowerBalance[t in T],
            sum(p[g, t] for g in G)  + sum(w[j,t] for j in J) - sum(demands_EL[d].D_el[t] for d in D_el) == 0 
        PowerFlowUp[l in L, t in T],
            lines[l].ptdf[[generators[g].ELnode for g in G]]'*([p[g, t] for g in G]) +
            lines[l].ptdf[[windgenerators[j].ELnode for j in J]]'*([w[j, t] for j in J]) -
            lines[l].ptdf[[demands_EL[d].ELnode for d in D_el]]'*([demands_EL[d].D_el[t] for d in D_el]) <= lines[l].fmax
        
        PowerFlowDown[l in L, t in T],
            lines[l].ptdf[[generators[g].ELnode for g in G]]'*([p[g, t] for g in G]) +
            lines[l].ptdf[[windgenerators[j].ELnode for j in J]]'*([w[j, t] for j in J]) -
            lines[l].ptdf[[demands_EL[d].ELnode for d in D_el]]'*([demands_EL[d].D_el[t] for d in D_el]) >= -lines[l].fmax 
    end
    
    if ES.print_model == 1
        print(m)
    end

    return m
end


# ####################### MISOCP #########################
function build_MISOCP_model(ES::IntegratedEnergySystem)

    M = ES.M # Big-M
    ########## Gas system ##########
    nodes, sources, demands, pipes, compressors,
        N, S, D, P, T, C, dt, C_cur_load_gas = get_gas_system_attributes(ES)

    ########## Power system ##########
    generators, windgenerators, demands_EL, lines,
        B, G, J, G, D_el, L = get_power_system_attributes(ES)
    GFPPs = [g for g in G if generators[g].type == "NGFPP"]
    nonGFPPs = [g for g in G if generators[g].type == "non-NGFPP"]

    #################### Initializing model ####################
    m = Model()
    set_solver(m, ES)

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
        Qd_gfpp[t in T, gg in GFPPs] # Gas consumption from GFPPs [kg/s]
        p[g in G, t in T] >= 0 # power dispatched by generators[MW]
        w[j in J, t in T] >= 0 # power dispatched by wind turbines [MW]
    end

    if isempty(compressors)==false # if there are compressors in the nework
        @variable(m, Q_c[t in T, c in C]) # Gas flow in each compressor [kg/s]
    end


    ################### Objective function ####################
    # Cost of gas supply (including gas used for GFPPs) and curtailment
    cost_gas = dt/3600 *(sum(Q_w[t,s]*sources[s].c1 + Q_w[t,s]^2*sources[s].c2 for s in S, t in T) +
        sum(Qd_cur_gas[t,d]*C_cur_load_gas for d in D, t in T)) #+
        # (isempty(compressors) ? 0.0 : sum(compressors[c].compr_cost*(
        #     pr[t,compressors[c].stop]-pr[t,compressors[c].start]) for c in C, t in T))
    # Cost of running non-GFPPs
    cost_el = dt/3600 *(sum(generators[g].c_el*p[g, t] for g in nonGFPPs, t in T))

    @objective(m, Min, cost_gas + cost_el) 

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
        sum(Q_nm_in⁻[t,p] for p in P if pipes[p].stop == n) -
        sum(generators[gg].eta*p[gg, t] for gg in GFPPs if generators[gg].NGnode == n))

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

    @constraints m begin
        GenerationLimitUp[g in G, t in T],
            p[g, t] <= generators[g].pmax
        # RampUp[g in G, t in T[2:end]], 
        #     p[g,t]-p[g, t-1] <= generators[g].pup/3600*dt
        # RampDown[g in G, t in T[2:end]],
        #     p[g,t]-p[g, t-1] >= -generators[g].pdown/3600*dt
        WindLimitUp[j in J, t in T],
            w[j, t] <= windgenerators[j].P_w[t]

        # Formulation with PTDF matrix
        PowerBalance[t in T],
            sum(p[g, t] for g in G)  + sum(w[j,t] for j in J) - sum(demands_EL[d].D_el[t] for d in D_el) == 0 
        PowerFlowUp[l in L, t in T],
            lines[l].ptdf[[generators[g].ELnode for g in G]]'*([p[g, t] for g in G]) +
            lines[l].ptdf[[windgenerators[j].ELnode for j in J]]'*([w[j, t] for j in J]) -
            lines[l].ptdf[[demands_EL[d].ELnode for d in D_el]]'*([demands_EL[d].D_el[t] for d in D_el]) <= lines[l].fmax
        
        PowerFlowDown[l in L, t in T],
            lines[l].ptdf[[generators[g].ELnode for g in G]]'*([p[g, t] for g in G]) +
            lines[l].ptdf[[windgenerators[j].ELnode for j in J]]'*([w[j, t] for j in J]) -
            lines[l].ptdf[[demands_EL[d].ELnode for d in D_el]]'*([demands_EL[d].D_el[t] for d in D_el]) >= -lines[l].fmax 
    end
    
    if ES.print_model == 1
        print(m)
    end

    return m
end


# ####################### Piecewise linear incremental #########################
function build_PWL_incremental_model(ES::IntegratedEnergySystem, no_V)

    ########## Gas system ##########
    nodes, sources, demands, pipes, compressors,
        N, S, D, P, T, C, dt, C_cur_load_gas = get_gas_system_attributes(ES)

    ########## Power system ##########
    generators, windgenerators, demands_EL, lines,
        B, G, J, G, D_el, L = get_power_system_attributes(ES)
    GFPPs = [g for g in G if generators[g].type == "NGFPP"]
    nonGFPPs = [g for g in G if generators[g].type == "non-NGFPP"]

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
    set_solver(m, ES)
    
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
        Qd_gfpp[t in T, gg in GFPPs] # Gas consumption from GFPPs [kg/s]
        p[g in G, t in T] >= 0 # power dispatched by generators[MW]
        w[j in J, t in T] >= 0 # power dispatched by wind turbines [MW]
    end
    
    if isempty(compressors)==false # if there are compressors in the nework
        @variable(m, Q_c[t in T, c in C]) # Gas flow in each compressor [kg/s]
    end


    ################### Objective function ####################
    # Cost of gas supply (including gas used for GFPPs) and curtailment
    cost_gas = dt/3600 *(sum(Q_w[t,s]*sources[s].c1 + Q_w[t,s]^2*sources[s].c2 for s in S, t in T) +
        sum(Qd_cur_gas[t,d]*C_cur_load_gas for d in D, t in T)) #+
        # (isempty(compressors) ? 0.0 : sum(compressors[c].compr_cost*(
        #     pr[t,compressors[c].stop]-pr[t,compressors[c].start]) for c in C, t in T))
    # Cost of running non-GFPPs
    cost_el = dt/3600 *(sum(generators[g].c_el*p[g, t] for g in nonGFPPs, t in T))

    @objective(m, Min, cost_gas + cost_el) 

    ################### Constraints ####################
    ############# Constraints Gas System #############
    @constraints m begin
        supply_limits[t in T, s in S],
            sources[s].Q_min <= Q_w[t,s] <= sources[s].Q_max
        pressure_limits[t in T, n in N],
            nodes[n].pr_min <= pr[t,n] <= nodes[n].pr_max
        avg_flow_def[t in T, p in P],
            Q_nm[t,p] == (Q_nm_in[t,p] + Q_nm_out[t,p])/2
        linepack_def[t in T, p in P],
            LP[t,p] == pipes[p].S*(pr[t,pipes[p].start]+pr[t,pipes[p].stop])/2
        linepack_intertemporal_initial[t = [1], p in P],
            #ϕ⁺[t,p] /2 / LP[t,p] * (
            Q_nm_in[t,p] - Q_nm_out[t,p] == 0
        linepack_intertemporal[t in T, p in P; t > 1],
            LP[t,p] == LP[t-1,p] +
                dt * (Q_nm_in[t,p] - Q_nm_out[t,p])
        linepack_end_condition[t = [T[end]], p in P],
            LP[t,p] >= LP[1,p] # pipes[p].LP0
        linearized_weymouth[t in T, p in P],
            Q_nm_sq_appr[t,p] == pipes[p].Knm^2 * (pr_sq_appr[t,pipes[p].start] - pr_sq_appr[t,pipes[p].stop])
        curtailment_limit[t in T, d in D],
            0 <= Qd_cur_gas[t,d] <= demands[d].Qd_gas[t]
    end

    if isempty(compressors) == false
        @constraints m begin
            compressor_flow[t in T, c in C],
                0 <= Q_c[t,c]
            compressor_limits_l[t in T, c in C],
                compressors[c].CR_min*pr[t,compressors[c].start] <= pr[t,compressors[c].stop] 
            compressor_limits_u[t in T, c in C],
                pr[t,compressors[c].stop] <= compressors[c].CR_max*pr[t,compressors[c].start]
        end
    end

    ###### Constraints linearization ######
    @constraints m begin
        Q_nm_sq_appr_def[t in T, p in P],
            Q_nm_sq_appr[t,p] ==  Q_set[t,p][1]*abs(Q_set[t,p][1]) +
                sum((Q_set[t,p][k+1]*abs(Q_set[t,p][k+1])-Q_set[t,p][k]*abs(Q_set[t,p][k]))*δ_Q[t,p,k] for k in K)
        Q_nm_appr_def[t in T, p in P],
            Q_nm[t,p] == Q_set[t,p][1] + sum((Q_set[t,p][k+1]-Q_set[t,p][k])*δ_Q[t,p,k] for k in K)
        δ_Q_limit[t in T, p in P, k in K[1:end-1]],
            δ_Q[t,p,k+1] <= y_Q[t,p,k]
        y_Q_limit[t in T, p in P, k in K[1:end-1]],
            y_Q[t,p,k] <= δ_Q[t,p,k] 
        pr_sq_appr_def[t in T, n in N],
            pr_sq_appr[t,n] ==  pr_set[t,n][1]*abs(pr_set[t,n][1]) +
                sum((pr_set[t,n][k+1]*abs(pr_set[t,n][k+1])-pr_set[t,n][k]*abs(pr_set[t,n][k]))*δ_pr[t,n,k] for k in K)
        pr_appr_def[t in T, n in N],
            pr[t,n] == pr_set[t,n][1] + sum((pr_set[t,n][k+1]-pr_set[t,n][k])*δ_pr[t,n,k] for k in K)
        δ_pr_limit[t in T, n in N, k in K[1:end-1]],
            δ_pr[t,n,k+1] <= y_pr[t,n,k]
        y_pr_limit[t in T, n in N, k in K[1:end-1]],
            y_pr[t,n,k] <= δ_pr[t,n,k] 
    end

    ################### Coupling Constraints ####################
    @expression(m, expr_node_balance[t in T, n in N],
        sum(Q_w[t,s] for s in S if sources[s].location == n) -
        sum(demands[d].Qd_gas[t] for d in D if demands[d].node_load == n) +
        sum(Qd_cur_gas[t,d] for d in D if demands[d].node_load == n) +
        sum(Q_nm_out[t,p] for p in P if pipes[p].stop == n) -
        sum(Q_nm_in[t,p] for p in P if pipes[p].start == n) -
        sum(generators[gg].eta*p[gg, t] for gg in GFPPs if generators[gg].NGnode == n))

    # Nodal balance equation with and without compressors
    if isempty(compressors) == true # if there are no compressors in the nework
        @constraint(m, node_balance[t in T, n in N], expr_node_balance[t,n] == 0)
    else 
        @constraint(m, node_balance[t in T, n in N],
            expr_node_balance[t,n] +
            sum(Q_c[t,c] for c in C if compressors[c].stop == n) -
            sum(Q_c[t,c] for c in C if compressors[c].start == n ) == 0)
    end

    ################### Constraints Power System ####################
    @constraints m begin
        GenerationLimitUp[g in G, t in T],
            p[g, t] <= generators[g].pmax
        # RampUp[g in G, t in T[2:end]], 
        #     p[g,t]-p[g, t-1] <= generators[g].pup/3600*dt
        # RampDown[g in G, t in T[2:end]],
        #     p[g,t]-p[g, t-1] >= -generators[g].pdown/3600*dt
        WindLimitUp[j in J, t in T],
            w[j, t] <= windgenerators[j].P_w[t]

        # Formulation with PTDF matrix
        PowerBalance[t in T],
            sum(p[g, t] for g in G)  + sum(w[j,t] for j in J) - sum(demands_EL[d].D_el[t] for d in D_el) == 0 
        PowerFlowUp[l in L, t in T],
            lines[l].ptdf[[generators[g].ELnode for g in G]]'*([p[g, t] for g in G]) +
            lines[l].ptdf[[windgenerators[j].ELnode for j in J]]'*([w[j, t] for j in J]) -
            lines[l].ptdf[[demands_EL[d].ELnode for d in D_el]]'*([demands_EL[d].D_el[t] for d in D_el]) <= lines[l].fmax
        
        PowerFlowDown[l in L, t in T],
            lines[l].ptdf[[generators[g].ELnode for g in G]]'*([p[g, t] for g in G]) +
            lines[l].ptdf[[windgenerators[j].ELnode for j in J]]'*([w[j, t] for j in J]) -
            lines[l].ptdf[[demands_EL[d].ELnode for d in D_el]]'*([demands_EL[d].D_el[t] for d in D_el]) >= -lines[l].fmax 
    end

    if ES.print_model == 1
        print(m)
    end

    return m

end


####################### MISOCP McCormick #########################
function build_MISOCP_McCormick_model(ES::IntegratedEnergySystem)

    M = ES.M # Big-M
    ########## Gas system ##########
    nodes, sources, demands, pipes, compressors,
        N, S, D, P, T, C, dt, C_cur_load_gas = get_gas_system_attributes(ES)

    ########## Power system ##########
    generators, windgenerators, demands_EL, lines,
        B, G, J, G, D_el, L = get_power_system_attributes(ES)
    GFPPs = [g for g in G if generators[g].type == "NGFPP"]
    nonGFPPs = [g for g in G if generators[g].type == "non-NGFPP"]

    m = Model()
    set_solver(m, ES)
 
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
        Qd_gfpp[t in T, gg in GFPPs] # Gas consumption from GFPPs [kg/s]
        p[g in G, t in T] >= 0 # power dispatched by generators[MW]
        w[j in J, t in T] >= 0 # power dispatched by wind turbines [MW]
    end

    if isempty(compressors)==false # if there are compressors in the nework
        @variable(m, Q_c[t in T, c in C]) # Gas flow in each compressor [kg/s]
    end

    ################### Objective function ####################
    # Cost of gas supply (including gas used for GFPPs) and curtailment
    cost_gas = dt/3600 *(sum(Q_w[t,s]*sources[s].c1 + Q_w[t,s]^2*sources[s].c2 for s in S, t in T) +
        sum(Qd_cur_gas[t,d]*C_cur_load_gas for d in D, t in T)) #+
        # (isempty(compressors) ? 0.0 : sum(compressors[c].compr_cost*(
        #     pr[t,compressors[c].stop]-pr[t,compressors[c].start]) for c in C, t in T))
    # Cost of running non-GFPPs
    cost_el = dt/3600 *(sum(generators[g].c_el*p[g, t] for g in nonGFPPs, t in T))

    @objective(m, Min, cost_gas + cost_el) 

    ################### Constraints ####################
    ############# Constraints Gas System #############
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

    if isempty(compressors) == false
        @constraints m begin
            compressor_flow[t in T, c in C],
                0 <= Q_c[t,c]
            compressor_limits_l[t in T, c in C],
                compressors[c].CR_min*pr[t,compressors[c].start] <= pr[t,compressors[c].stop] 
            compressor_limits_u[t in T, c in C],
                pr[t,compressors[c].stop] <= compressors[c].CR_max*pr[t,compressors[c].start]
        end
    end

    ################### Coupling Constraints ####################
    @expression(m, expr_node_balance[t in T, n in N],
        sum(Q_w[t,s] for s in S if sources[s].location == n) -
        sum(demands[d].Qd_gas[t] for d in D if demands[d].node_load == n) +
        sum(Qd_cur_gas[t,d] for d in D if demands[d].node_load == n) +
        sum(Q_nm_out⁺[t,p] for p in P if pipes[p].stop == n) -
        sum(Q_nm_in⁺[t,p] for p in P if pipes[p].start == n) +
        sum(Q_nm_out⁻[t,p] for p in P if pipes[p].start == n) -
        sum(Q_nm_in⁻[t,p] for p in P if pipes[p].stop == n) -
        sum(generators[gg].eta*p[gg, t] for gg in GFPPs if generators[gg].NGnode == n))

    # Nodal balance equation with and without compressors
    if isempty(compressors) == true # if there are no compressors in the nework
        @constraint(m, node_balance[t in T, n in N], expr_node_balance[t,n] == 0)
    else 
        @constraint(m, node_balance[t in T, n in N],
            expr_node_balance[t,n] +
            sum(Q_c[t,c] for c in C if compressors[c].stop == n) -
            sum(Q_c[t,c] for c in C if compressors[c].start == n ) == 0)
    end

    ################### Constraints Power System ####################
    @constraints m begin
        GenerationLimitUp[g in G, t in T],
            p[g, t] <= generators[g].pmax
        # RampUp[g in G, t in T[2:end]], 
        #     p[g,t]-p[g, t-1] <= generators[g].pup/3600*dt
        # RampDown[g in G, t in T[2:end]],
        #     p[g,t]-p[g, t-1] >= -generators[g].pdown/3600*dt
        WindLimitUp[j in J, t in T],
            w[j, t] <= windgenerators[j].P_w[t]

        # Formulation with PTDF matrix
        PowerBalance[t in T],
            sum(p[g, t] for g in G)  + sum(w[j,t] for j in J) - sum(demands_EL[d].D_el[t] for d in D_el) == 0 
        PowerFlowUp[l in L, t in T],
            lines[l].ptdf[[generators[g].ELnode for g in G]]'*([p[g, t] for g in G]) +
            lines[l].ptdf[[windgenerators[j].ELnode for j in J]]'*([w[j, t] for j in J]) -
            lines[l].ptdf[[demands_EL[d].ELnode for d in D_el]]'*([demands_EL[d].D_el[t] for d in D_el]) <= lines[l].fmax
        
        PowerFlowDown[l in L, t in T],
            lines[l].ptdf[[generators[g].ELnode for g in G]]'*([p[g, t] for g in G]) +
            lines[l].ptdf[[windgenerators[j].ELnode for j in J]]'*([w[j, t] for j in J]) -
            lines[l].ptdf[[demands_EL[d].ELnode for d in D_el]]'*([demands_EL[d].D_el[t] for d in D_el]) >= -lines[l].fmax 
    end

    if ES.print_model == 1
        print(m)
    end

    return m 
end