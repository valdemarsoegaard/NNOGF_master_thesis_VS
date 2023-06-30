include("energy_system.jl")
include("PTDF_matrix.jl")

##================= Gas nodes data ===============
function create_Nodes(path_input_data)
    N_df = CSV.read(joinpath(path_input_data, "gas_nodes.csv"), DataFrame, delim=",")
    nodes = Dict()
    for row in 1:nrow(N_df)
        id, pr_max, pr_min, node_type =
            N_df[row, [:Node_No, :Pmax_MPa, :Pmin_MPa, :Node_Type]]
        nodes[id] =  Node(id, pr_max, pr_min, node_type) 
    end
    return nodes
end


##================= Gas sources data ===============
function create_GasSources(path_input_data)
    S_df = CSV.read(joinpath(path_input_data, "gas_supply.csv"), DataFrame, delim=",")
    sources = Dict()
    for row in 1:nrow(S_df)
        id, location, Q_max, Q_min, c1, c2 =
            S_df[row, [:Supply_No, :Node, :Smax_kg_s, :Smin_kg_s,:C1_per_kgh, :C2_per_kgh2]]
        sources[id] =  GasSource(id, location, Q_max, Q_min, c1, c2) 
    end
    return sources
end


##================= Gas load data===============
function create_GasDemands(path_input_data, dt, N_ts_h)
    D_df = CSV.read(joinpath(path_input_data, "gas_load.csv"), DataFrame, delim=",")
    demands = Dict()
    GP_df= CSV.read(joinpath(path_input_data, "gas_profile.csv"), DataFrame, delim=",")
    Num_data= nrow(GP_df) # Number of input load data
    if Num_data==24 #Uniform gas model: constant load every hour
        unif=1;
    elseif Num_data==480 # Data provided every 3 min (60*24/3) and average based on dt chosen by user)
        unif=0;
        dt_data=3*60; # 180 seconds in 3 minutes
    elseif Num_data==288 # Data provided every 5 min (60*24/5) and average based on dt chosen by user)
        unif=0;
        dt_data=5*60; # 300 seconds in 5 minutes
    elseif Num_data==96 # Data provided every 15 min (60*24/15) and average based on dt chosen by user)
        unif=0;
        dt_data=15*60; # 300 seconds in 5 minutes
    else
        throw(ArgumentError("Input gas load data must have either 3, 5, 15, or 60 minute resolution."))
    end

    for row in 1:nrow(D_df)
        id, node_load, load_nom, type_profile = 
            D_df[row, [:Load_No,:Node,:Load_kg_s,:Profile]]
        load_profile = GP_df[:,type_profile]
        if  unif==1
            demands[id] =  GasLoadUniform(
                id, node_load, load_nom, type_profile, load_profile, N_ts_h) 
        else
            demands[id] =  GasLoadAverage(
                id, node_load, load_nom, type_profile, load_profile, dt, dt_data) 
        end
    end
    return demands
end


##================= Gas pipeline data ===============
function create_Pipes(path_input_data, nodes, cc)
    P_df = CSV.read(joinpath(path_input_data, "gas_pipes.csv"), DataFrame, delim=",")
    pipes = Dict()
    for row in 1:nrow(P_df)
        id, start, stop, fr, DD, LL =
            P_df[row, [:Pipe_No, :From_Node, :To_Node, :friction, :Diameter_m,:Length_m]]
         pipes[id] =  Pipeline(id, start, stop, fr, DD, LL, cc, nodes) 
    end
    return pipes
end


##================= Gas compressor data ===============
function create_Compressors(path_input_data, nodes)
    C_df = CSV.read(joinpath(path_input_data, "gas_compressors.csv"), DataFrame, delim=",")
    compressors = Dict()
    for row in 1:nrow(C_df)
        id, from_node, to_node, CR_max, CR_min, compr_cost =
            C_df[row, [:Compressor_No,:From_Node,:To_Node,:CR_Max,:CR_Min,:Compression_cost]]
        compressors[id] =  Compressor(
            id, from_node, to_node, CR_max, CR_min, compr_cost, nodes) 
    end
    return compressors
end


##================= Electricity bus data ===============
function create_Buses(path_input_data)
    # CSV.read(joinpath(dirname(@__DIR__), "inputs", "electricity",case_study_el,"buses_EL.csv"), DataFrame, delim=",")
    B_df = CSV.read(joinpath(path_input_data, "buses_EL.csv"), DataFrame, delim=",")
    buses = Dict()
    for row in 1:nrow(B_df)
        id, slack  = B_df[row, [:Bus_No, :Slack]]
        buses[id] =  Bus(id, slack) 
    end
    return buses
end


##================= Electricity DispatchableGenerators data ===============
function create_DispatchableGenerators(path_input_data)
    G_df = CSV.read(joinpath(path_input_data, "dispatchablegenerators.csv"), DataFrame, delim=",")
    generators = Dict()
    for row in 1:nrow(G_df)
        id, ELnode, pmin, pmax, pdown, pup, type, NGnode, eta, c_el  = G_df[row, [:Gen_num,:EL_node,:Pmin_MW,:Pmax_MW,:P_down_MW_h,:P_up_MW_h,:Type,:NG_node,:Conversion_kg_sMW,:Cost_per_MWh]]
        generators[id] =  DispatchableGenerator(id, ELnode, pmin, pmax, pdown, pup, type, NGnode, eta, c_el ) 
    end
    return generators
end


##================= Electricity WindGenerators data ===============
function create_WindGenerators(path_input_data, dt)
    #CSV.read(joinpath(dirname(@__DIR__), "inputs","electricity",case_study_el,"windgenerators.csv"), DataFrame, delim=",")
    #CSV.read(joinpath(dirname(@__DIR__), "inputs","electricity",case_study_el,"wind_profile.csv"), DataFrame, delim=",")
    W_df = CSV.read(joinpath(path_input_data, "windgenerators.csv"), DataFrame, delim=",")
    windgenerators = Dict()
    WP_df= CSV.read(joinpath(path_input_data, "wind_profile.csv"), DataFrame, delim=",")
    Num_data_w= nrow(WP_df)
    T_data=24
    dt_data=T_data*3600/Num_data_w; # input profile provided every dt_data seconds
        for row in 1:nrow(W_df)
            id, ELnode, p_max, type_profile = W_df[row, [:Wind_num,:EL_node,:Pmax_MW,:profile_type]]
            wind_profile = WP_df[:,type_profile]
            windgenerators[id] =  WindGenerator(id, ELnode, p_max, type_profile, wind_profile, dt, dt_data) 
        end
    return windgenerators
end


##================= Electricity demand data ===============
#CSV.read(joinpath(dirname(@__DIR__), "inputs","electricity",case_study_el,"electricity_load.csv"), DataFrame, delim=",")
#CSV.read(joinpath(dirname(@__DIR__), "inputs","electricity",case_study_el,"electricity_profile.csv"), DataFrame, delim=",")
function create_PowerDemands(path_input_data, dt)
    D_el_df = CSV.read(joinpath(path_input_data, "electricity_load.csv"), DataFrame, delim=",")
    demand_EL = Dict()
    EP_df= CSV.read(joinpath(path_input_data, "electricity_profile.csv"), DataFrame, delim=",")
    Num_data_el= nrow(EP_df) # Number of input load data
    T_data=24
    dt_data=T_data*3600/Num_data_el; # input profile provided every dt_data seconds
    for row in 1:nrow(D_el_df)
        id, ELnode, load_nom, type_profile = D_el_df[row, [:Load_No,:EL_Node,:Load_MW,:Profile]]
        load_profile = EP_df[:,type_profile]
        demand_EL[id] =  PowerDemand(id, ELnode, load_nom, type_profile, load_profile, dt, dt_data) 
    end
    return demand_EL
end


##================= Electricity Line data ===============
#dirname(@__DIR__), "inputs","electricity",case_study_el,
function create_Lines(path_input_data)
    L_df = CSV.read(joinpath(path_input_data, "lines.csv"), DataFrame, delim=",")
    Ψ = create_PTDF_matrix(L_df)
    Ψ = round.(Ψ, digits=4)

    lines= Dict()
    for row in 1:nrow(L_df)
        id, start, stop, x, fmax = L_df[row, [:Line_num,:Start,:Stop,:X_pu,:Capacity_MW]]
        lines[id] =  Line(id, start, stop, x, fmax, Ψ[row,:]) 
    end
    return lines
end


function intialize_data!(GS :: GasSystem)
    """
    Loads and returns input data from CSV files for a given case study.
    """
    ##================= Global parameters ===============
    path_input_data_gas = GS.path_input_data_gas

    ##================= Instantiate gas network data ===============
    nodes = create_Nodes(path_input_data_gas)
    pipes = create_Pipes(path_input_data_gas, nodes, GS.cc)
    GS.sources = create_GasSources(path_input_data_gas) 
    GS.demands = create_GasDemands(path_input_data_gas, GS.dt, GS.N_ts_h)
    GS.compressors = create_Compressors(path_input_data_gas, nodes)

    ##================= Run space discretization gas network ===============
    if GS.space_disc == 1
        GS.nodes, GS.pipes = space_discretization(
            nodes, pipes, GS.segment)
    else
        GS.nodes, GS.pipes = nodes, pipes
    end

    ##================= Create ordered sets ===============
    GS.N = sort(collect(keys(GS.nodes)))
    GS.S = sort(collect(keys(GS.sources)))
    GS.P = sort(collect(keys(GS.pipes)))
    GS.D = sort(collect(keys(GS.demands)))
    GS.C = sort(collect(keys(GS.compressors)))

end


function intialize_data!(IES :: IntegratedEnergySystem)
    """
    Loads and returns input data from CSV files for a given case study.
    """
    ##================= Global parameters ===============
    path_input_data_gas = IES.path_input_data_gas
    path_input_data_power = IES.path_input_data_power


    ##================= Instantiate gas network data ===============
    nodes = create_Nodes(path_input_data_gas)
    pipes = create_Pipes(path_input_data_gas, nodes, IES.cc)
    IES.sources = create_GasSources(path_input_data_gas) 
    IES.demands = create_GasDemands(path_input_data_gas, IES.dt, IES.N_ts_h)
    IES.compressors = create_Compressors(path_input_data_gas, nodes)

    ##================= Run space discretization gas network ===============
    if IES.space_disc == 1
        IES.nodes, IES.pipes = space_discretization(
            nodes, pipes, IES.segment)
    else
        IES.nodes, IES.pipes = nodes, pipes
    end

    ##================= Instantiate power network data ===============
    IES.buses = create_Buses(path_input_data_power)
    IES.lines = create_Lines(path_input_data_power)
    IES.generators = create_DispatchableGenerators(path_input_data_power)
    IES.windgenerators = create_WindGenerators(path_input_data_power, IES.dt)
    IES.demands_EL = create_PowerDemands(path_input_data_power, IES.dt)

    ##================= Create ordered sets ===============
    ##### Gas #####
    IES.N = sort(collect(keys(IES.nodes)))
    IES.S = sort(collect(keys(IES.sources)))
    IES.P = sort(collect(keys(IES.pipes)))
    IES.D = sort(collect(keys(IES.demands)))
    IES.C = sort(collect(keys(IES.compressors)))

    ##### Power #####
    IES.B = sort(collect(keys(IES.buses)))
    IES.L = sort(collect(keys(IES.lines)))
    IES.G = sort(collect(keys(IES.generators)))
    IES.J = sort(collect(keys(IES.windgenerators)))
    IES.D_el = sort(collect(keys(IES.demands_EL)))

end


function PDE_type_parameters(PDE_type::Int)
    """
    Chooses model parameters for the PDEs of the gas flow model.
    """
    if PDE_type==1 ut1=0; ut2=0; #steady-state model
    elseif PDE_type==2 ut1=1; ut2=0; #quasi-dynamic model
    elseif PDE_type==3 ut1=1; ut2=1; #transient model
    end

    return ut1, ut2
end

function read_NN_parameters(NN_data::Dict{String,Any}) #read a .json file with the weights and biases of the neural network
    """
    Reads the weights and biases of the neural network from a dictionary created by a .json file.
    """

    unsqueeze(a::Array) = begin  #functions to unsqueeze a vector to a matrix. This is necesarry due to the .json format
        nsize = foldl(append!, size(a); init=[1])
        return reshape(a, Tuple(nsize))
    end
    stack(vs::Vector{T}) where T<:Vector = begin
        return reduce(vcat, map(unsqueeze, map(stack, vs)))
    end
    stack(vs::Vector{T}) where T<:Real = vs

    K=sum([occursin("weight", k) for k in keys(NN_data)]) #count how many layers are in the neural network


    #Save weights in a dictionary for each layer
    W = Dict{Int64,Any}()
    for i in 1:K
    # W[i] = stack([convert(Array{Float64,1}, NN_data["layer$i.weight"][k]) for k in 1:size(NN_data["layer$i.weight"])[1]])
    W[i] = stack([convert(Array{Float64,1}, NN_data["model.layer$i.weight"][k]) for k in 1:size(NN_data["model.layer$i.weight"])[1]])
    end
   

    #Save biases in a dictionary for each layer
    B = Dict{Int64,Any}()
    for i in 1:K
        if typeof(NN_data["model.layer$i.bias"]) == Array{Any,1}
            B[i] = convert(Array{Float64,1}, NN_data["model.layer$i.bias"]) 
        else
            B[i] = NN_data["model.layer$i.bias"]
        end
    end

    #load leaky relu parameter (if any)
    if "α" in keys(NN_data)
        α = NN_data["α"]
    else
        α = 0
    end 
    
    nn_data_keys = keys(NN_data)
    #load bounds of the NN
    y1_bounds = stack([ convert(Array{Float64,1}, NN_data["y1_hat"][k]) for k in 1:size(NN_data["y1_hat"])[1]] )
    y_hat_min = y1_bounds[1,:]
    y_hat_max = y1_bounds[2,:]
    if "y2_hat" in nn_data_keys
        y2_bounds = stack([ convert(Array{Float64,1}, NN_data["y2_hat"][k]) for k in 1:size(NN_data["y2_hat"])[1]] )
        y_hat_min = [y1_bounds[1,:] y2_bounds[1,:]]
        y_hat_max = [y1_bounds[2,:] y2_bounds[2,:]]
    end
    if "y3_hat" in nn_data_keys
        y3_bounds = stack([ convert(Array{Float64,1}, NN_data["y3_hat"][k]) for k in 1:size(NN_data["y3_hat"])[1]] )
        y_hat_min = [y1_bounds[1,:] y2_bounds[1,:] y3_bounds[1,:]]
        y_hat_max = [y1_bounds[2,:] y2_bounds[2,:] y3_bounds[2,:]]
    end
   
    

    #load data scaler parameters
    mn = NN_data["mean"]
    sd = NN_data["std"]


    return W, B, 1:K-1, α, y_hat_min, y_hat_max, mn, sd
end