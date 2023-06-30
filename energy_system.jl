abstract type EnergySystem end

mutable struct GasSystem <: EnergySystem
    cc #:: Int64 # Speed of sound 350 m/s
    path_input_data_gas :: String # input folder of scenario data
    case_study :: String # Case study
    nodes #:: Dict{Int64, Node}() # Node container
    sources #:: Dict{Int64, Source}() # Source container
    demands #:: Dict{Int64, Any}() # Demand container -- can be either type Gasload_uniform or Gasload_average
    pipes #:: Dict{Int64, Pipeline}() # Pipe container
    compressors #:: Dict{Int64, Compressor}() # Compressor container
    N # Set of nodes
    S # Set of sources
    P # Set of pipes
    D # Set of demands
    C # Set of compressors
    space_disc :: Int64 # Control parameter if space discretization should be performed
    segment :: Int64 # Length of each subpipeline segment
    timehorizon :: Int64 # Time horizon (hours)
    print_model :: Int64 # Controls if model is printed
    T # Set of timesteps
    no_V # Number of set points for linearization (PWL)
    dt :: Int64 # Time interval for discretization (seconds)
    N_dt :: Int64 # Number of timesteps in the time horizon#
    N_ts_h :: Int64 # Number of timesteps per hour
    # ut1 :: Int64 # Control parameter for PDE model
    # ut2 :: Int64 # Control parameter for PDE model
    PDE_type :: Int64 # Control parameter for PDE model
    model_type :: String # Model type
    C_cur_load_gas :: Int64 # Value of lost load (\$/kgh)
    M # Parameter for Big-M method
    obj_val # Objective function value of model
    solve_time # Solve time of model instance
    Q_w # Optimal gas production
    pr # Optimal nodal pressure
    Q_nm # Optimal gas flow
    Q_nm_in # Optimal gas flow in
    Q_nm_out # Optimal gas flow Q_nm_out
    pr_avg # Optimal average pressure in pipelines
    Qd_cur_gas # Optimal curtailment
    Q_c # Optimal gas flow through compressors
    LP # Optimal linepack mass
    ϕ⁺ # Optimal value of auxilary variable related to pressures
    ϕ⁻ # Optimal value of auxilary variable related to pressures
    y_nm # Optimal flow directions
    Q_nm⁺ # Optimal average flow in pipeline nm in postive direction
    Q_nm_in⁺ # Optimal inflow in pipeline nm in postive direction
    Q_nm_out⁺ # Optimal outflow in pipeline nm in postive direction
    Q_nm⁻ # Optimal average flow in pipeline nm in postive direction
    Q_nm_in⁻ # Optimal inflow in pipeline nm in postive direction
    Q_nm_out⁻ # Optimal outflow in pipeline nm in postive direction
    pr_restored # restored pressure values from original PDE problem
    obj_val_restored # restored objective function value from original PDE problem
    config_dict_solver # Tuple, contains (1) solver and (2) solver attribute dictionary
    model # JuMP model
    out_dic # Output decision variables
    obj_type :: String #"linear" or "quadratic"
    callback :: Bool
    # Adapted from https://discourse.julialang.org/t/default-value-of-some-fields-in-a-mutable-struct/33408/20
    # Authors: GunnarFarneback, davidbp
    function GasSystem(timehorizon, dt, case_study; kwargs...)
        GS = new()
        GS.cc = 350
        GS.N_dt = Int64(timehorizon*3600/dt)
        GS.N_ts_h = Int64(3600/dt)
        GS.T = sort(collect(1:GS.N_dt))
        GS.path_input_data_gas = joinpath(dirname(@__DIR__), "inputs", "gas_only", case_study)
        for (key, value) in kwargs
            # field_type_key = typeof(getfield(GS, key))
            setfield!(GS, key, value)
            # setfield!(GS, key, convert(field_type_key, value))
        end
        return GS
    end
end


function get_gas_system_attributes(ES::EnergySystem)

    return ES.nodes, ES.sources, ES.demands, ES.pipes, ES.compressors,
        ES.N, ES.S, ES.D, ES.P, ES.T, ES.C, ES.dt, ES.C_cur_load_gas
end



mutable struct IntegratedEnergySystem <: EnergySystem
    ##### Gas system attributes #####
    cc :: Int64 # Speed of sound 350 m/s
    path_input_data_gas :: String # input folder of gas scenario data
    case_study :: String # Case study
    nodes #:: Dict{Int64, Node}() # Node container
    sources #:: Dict{Int64, Source}() # Source container
    demands #:: Dict{Int64, Any}() # Demand container -- can be either type Gasload_uniform or Gasload_average
    pipes #:: Dict{Int64, Pipeline}() # Pipe container
    compressors #:: Dict{Int64, Compressor}() # Compressor container
    N # Set of nodes
    S # Set of sources
    P # Set of pipes
    D # Set of demands
    C # Set of compressors
    space_disc :: Int64 # Control parameter if space discretization should be performed
    segment :: Int64 # Length of each subpipeline segment
    timehorizon :: Int64 # Time horizon (hours)
    print_model :: Int64 # Controls if model is printed
    T # Set of timesteps
    no_V # Number of set points for linearization (PWL)
    dt :: Int64 # Time interval for discretization (seconds)
    N_dt :: Int64 # Number of timesteps in the time horizon#
    N_ts_h :: Int64 # Number of timesteps per hour
    # ut1 :: Int64 # Control parameter for PDE model
    # ut2 :: Int64 # Control parameter for PDE model
    PDE_type :: Int64 # Control parameter for PDE model
    model_type :: String # Model type
    C_cur_load_gas :: Int64 # Cost of lost load (\$/kgh)
    M # Parameter for Big-M method
    obj_val # Objective function value of model
    solve_time # Solve time of model instance
    Q_w # Optimal gas production
    pr # Optimal nodal pressure
    Q_nm # Optimal gas flow
    Q_nm_in # Optimal gas flow in
    Q_nm_out # Optimal gas flow Q_nm_out
    pr_avg # Optimal average pressure in pipelines
    Qd_cur_gas # Optimal curtailment
    Q_c # Optimal gas flow through compressors
    LP # Optimal linepack mass
    ϕ⁺ # Optimal value of auxilary variable related to pressures
    ϕ⁻ # Optimal value of auxilary variable related to pressures
    y_nm # Optimal flow directions
    Q_nm⁺ # Optimal average flow in pipeline nm in postive direction
    Q_nm_in⁺ # Optimal inflow in pipeline nm in postive direction
    Q_nm_out⁺ # Optimal outflow in pipeline nm in postive direction
    Q_nm⁻ # Optimal average flow in pipeline nm in postive direction
    Q_nm_in⁻ # Optimal inflow in pipeline nm in postive direction
    Q_nm_out⁻ # Optimal outflow in pipeline nm in postive direction
    pr_restored # restored pressure values from original PDE problem
    obj_val_restored # restored objective function value from original PDE problem

    ##### Power system attributes #####
    path_input_data_power :: String # input folder of gas scenario data
    buses # Dict containing power bus objects
    lines # Dict containing power line objects
    generators # Dict containing dispatchable generator objects
    windgenerators # Dict containing wind power objects
    demands_EL # Dict containing power demand objects
    B # List of buses
    L # List of lines
    G # List of dispatchable generators
    J # List of wind power producers
    D_el # List of power demands
    Qd_gfpp # Gas consumption of GFPPs at optimum
    p # Generator dispatch at optimum
    w # Wind generator dispatch at optimum
    config_dict_solver # Tuple, contains (1) solver and (2) solver attribute dictionary

    # Adapted from https://discourse.julialang.org/t/default-value-of-some-fields-in-a-mutable-struct/33408/20
    # Authors: GunnarFarneback, davidbp
    function IntegratedEnergySystem(timehorizon, dt, case_study; kwargs...)
        IES = new()
        IES.cc = 350
        IES.N_dt = Int64(timehorizon*3600/dt)
        IES.N_ts_h = Int64(3600/dt)
        IES.T = sort(collect(1:IES.N_dt))
        IES.path_input_data_gas = joinpath(dirname(@__DIR__), "inputs", "integrated_power_and_gas", case_study, "gas")
        IES.path_input_data_power = joinpath(dirname(@__DIR__), "inputs", "integrated_power_and_gas", case_study, "power")
        for (key, value) in kwargs
            # field_type_key = typeof(getfield(ES, key))
            setfield!(IES, key, value)
            # setfield!(ES, key, convert(field_type_key, value))
        end
        return IES
    end
end

function get_power_system_attributes(ES::IntegratedEnergySystem)

    return ES.generators, ES.windgenerators, ES.demands_EL, ES.lines,
        ES.B, ES.G, ES.J, ES.G, ES.D_el, ES.L
end


function intialize_energy_system(config_dict, system)

    if system == "gas_only"
        ES = GasSystem(
            config_dict[:timehorizon], 
            config_dict[:dt],
            config_dict[:case_study];
            config_dict...)

    elseif system == "integrated_power_and_gas"
        ES = IntegratedEnergySystem(
            config_dict[:timehorizon], 
            config_dict[:dt],
            config_dict[:case_study];
            config_dict...)
    else 
        throw(ArgumentError("System is set to $system.
            Must be either 'integrated_power_and_gas' or 'gas_only'."))
    end
    return ES
end