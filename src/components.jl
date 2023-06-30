############### Gas network structs ###############
mutable struct GasSource
    id # id
    location # Gas source location
    Q_max # Maximum source kg/s
    Q_min # Minimum source kg/s
    c1 # Linear coefficient, cost per kgh
    c2 # Quadratic coefficient, cost per kgh^2
end

mutable struct Node
    id # id
    pr_max # Maximum nodal pressure [MPa]
    pr_min # Minimum nodal pressure [MPa]
    node_type # slack nodes =1
end

mutable struct Pipeline
    id # id
    start # Pipeline start
    stop # Pipeline end
    fr # friction coefficient
    cc # # Speed of sound (m/s)
    DD # Diameter (m)
    LL # Lenght (m)
    AA # pipeline cross sectional area (m2)
    S # Linepack parameter
    LP0 # initial line pack mass
    Knm # Weymouth coefficient (kg/s)/MPa
    Q_nm_min # Min flow rate in pipelines (kg/s)
    Q_nm_max # Max flow rate in pipelines (kg/s)
    Q_nm_min_k # Array for storing updated minimum bounds in bound tightening method (kg/s)
    Q_nm_max_k # Array for storing updated maximum bounds in bound tightening method (kg/s)
    ϕ⁺_min :: Float64 # Bounds for auxiliary variable related to pressures (MPa)
    ϕ⁺_max :: Float64 # Bounds for auxiliary variable related to pressures (MPa)
    ϕ⁻_min :: Float64 # Bounds for auxiliary variable related to pressures (MPa)
    ϕ⁻_max :: Float64 # Bounds for auxiliary variable related to pressures (MPa)
    ϕ⁺_min_k #:: Vector{Vector{Float64}} # Array for storing updated auxiliary pressure bounds (MPA)
    ϕ⁺_max_k #:: Vector{Vector{Float64}} # Array for storing updated auxiliary pressure bounds (MPA)
    ϕ⁻_min_k #:: Vector{Vector{Float64}} # Array for storing updated auxiliary pressure bounds (MPA)
    ϕ⁻_max_k #:: Vector{Vector{Float64}} # Array for storing updated auxiliary pressure bounds (MPA)
    V_m
    V_p
    V_f
    function Pipeline(id, start, stop, fr, DD, LL, cc, nodes)
        AA = pi*(DD/2)^2
        Knm = sqrt(DD*AA^2/(fr*cc^2*LL))*10^6
        Q_nm_min = -Knm*sqrt(nodes[stop].pr_max^2-nodes[start].pr_min^2)
        Q_nm_max = Knm*sqrt(nodes[start].pr_max^2-nodes[stop].pr_min^2)
        ϕ⁺_min = nodes[start].pr_min + nodes[stop].pr_min
        ϕ⁺_max = nodes[start].pr_max + nodes[stop].pr_max
        ϕ⁻_min = nodes[start].pr_min - nodes[stop].pr_max
        ϕ⁻_max = nodes[start].pr_max - nodes[stop].pr_min
        V_m = cc^2/AA*10^-6
        V_p = AA*10^6
        V_f = fr*cc^2/(2*DD*AA)*10^-6 ### MPa!!!
        S = LL/V_m
        LP0 = LL/V_m*1 # 1 is arbitrary for average pressure
        # ϕ⁺_min_k = [[ϕ⁺_min for t in T]]
        # ϕ⁺_max_k = [[ϕ⁺_max for t in T]]
        # ϕ⁻_min_k = [[ϕ⁻_min for t in T]]
        # ϕ⁻_max_k = [[ϕ⁻_max for t in T]]
        return new(
            id, start, stop, fr, cc, DD, LL, AA, S, LP0, Knm, Q_nm_min, Q_nm_max, nothing, nothing,
            ϕ⁺_min, ϕ⁺_max, ϕ⁻_min, ϕ⁻_max, nothing, nothing, nothing, nothing, V_m, V_p, V_f)
    end
end

mutable struct Compressor
    id # id
    start # Start node
    stop # End node
    CR_max # Maximum compression ratio
    CR_min # Minimum compression ratio
    compr_cost # Compression cost (cost per unit of pressure increase)
    ϕ⁺_min_c # Bounds for auxiliary variable related to pressures (MPa)
    ϕ⁺_max_c # Bounds for auxiliary variable related to pressures (MPa)
    ϕ⁻_min_c # Bounds for auxiliary variable related to pressures (MPa)
    ϕ⁻_max_c # Bounds for auxiliary variable related to pressures (MPa)
    function Compressor(id, start, stop, CR_max, CR_min, compr_cost, nodes)
        ϕ⁺_min_c = nodes[start].pr_min + nodes[stop].pr_min
        ϕ⁺_max_c = nodes[start].pr_max + nodes[stop].pr_max
        ϕ⁻_min_c = nodes[start].pr_min - nodes[stop].pr_max
        ϕ⁻_max_c = nodes[start].pr_max - nodes[stop].pr_min
        return new(
            id, start, stop, CR_max, CR_min, compr_cost, 
            ϕ⁺_min_c, ϕ⁺_max_c, ϕ⁻_min_c, ϕ⁻_max_c)
    end
end

mutable struct GasLoadUniform
    id # id
    node_load # Node where the load is connected
    load_nom # Nominal load at each node kg/s
    type_profile # Type of gas load profile
    load_profile # Gas load profile every hour (-)
    Qd_gas # Gas load in kg/s at each timestep
    function GasLoadUniform(id,node_load, load_nom, type_profile, load_profile, N_ts_h)
        Qd_gas = repeat(round.(load_nom.*load_profile, digits=4), inner=N_ts_h)
        return new(id,node_load, load_nom, type_profile, load_profile, Qd_gas)
    end
end

mutable struct GasLoadAverage
    id # id
    node_load # Node where the load is connected
    load_nom # Nominal load at each node kg/s
    type_profile # Type of gas load profile
    load_profile # Gas load profile every 3 minutes (-)
    Qd_gas # Gas load in kg/s at each timestep
    function GasLoadAverage(id,node_load, load_nom, type_profile, load_profile, dt, dt_data)
        N_ave=Int64(dt/dt_data) # number of input data provided in timestep
        Qd_gas_dt = round.(load_nom.*load_profile, digits=4)
        Qd_gas = mean(reshape(Qd_gas_dt, N_ave, :), dims=1)
        return new(id,node_load, load_nom, type_profile, load_profile, Qd_gas)
    end
end


############### Power network structs ###############
## Components EL network
mutable struct Bus
    id # id
    slack # 1:slack bus, 0:otherwise
end

mutable struct DispatchableGenerator
    id # id
    ELnode # Generator location in the power system
    pmin # Minimum power [MW]
    pmax # Maximum power [MW]
    pdown # Maximmum ramp down rate [MW/h]
    pup # Maximmum ramp up rate [MW/h]
    type # Type of generator ("NGFPP" or "non-NGFPP")
    NGnode # Generator location in the gas system (for NGFPPs)
    eta # Power conversion factor (kcf/MWh) ## CHANGE UNIT OF MEASURE!!
    c_el # Linear coefficient, cost per MWh (for non-NGFPPs)
end

mutable struct WindGenerator
    id # id
    ELnode # Generator location in the power system
    pmax # Maximum power [MW]
    profile_type # Type of wind profile
    wind_profile # Daily wind profile in p.u.
    P_w # Wind generation in each timestep [MW]
    function WindGenerator(id,ELnode, pmax, profile_type, wind_profile, dt, dt_data)
        N_ave=Int64(dt/dt_data) # number of input data provided in timestep
        P_w_dt = round.(pmax.*wind_profile, digits=4)
        P_w = mean(reshape(P_w_dt, N_ave, :), dims=1)
        return new(id,ELnode, pmax, profile_type, wind_profile, P_w)
    end
end

mutable struct PowerDemand
    id # id
    ELnode # Node where the load is connected
    load_nom # Nominal electrical load [MW]
    type_profile # Type of gas load profile
    load_profile # Electrical load profile
    D_el # Electricity demand at each timestep
    function PowerDemand(id,ELnode, load_nom, type_profile, load_profile, dt, dt_data)
        N_ave=Int64(dt/dt_data) # number of input data provided in timestep
        D_el_dt = round.(load_nom.*load_profile, digits=4)
        D_el = mean(reshape(D_el_dt, N_ave, :), dims=1)
        return new(id,ELnode, load_nom, type_profile, load_profile, D_el)
    end
end

mutable struct Line
    id # id
    start # Line start
    stop # Line ends
    x # line reactance in p.u.
    fmax # Maximum capacity [MW]
    ptdf # For each line, row of the PTDF matrix
end
