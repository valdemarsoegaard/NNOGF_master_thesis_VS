using Pkg
Pkg.activate(".")
using JuMP, Ipopt, Gurobi, CSV, DataFrames, Statistics, Plots, LinearAlgebra, Dates, Tables
using Random
include("components.jl")
include("discretization.jl")
include("opt_models_gas.jl")
include("opt_models_integrated_power_and_gas.jl")
include("input_data.jl")
include("algorithms.jl")
include("energy_system.jl")
include("PTDF_matrix.jl")
include("model_functions.jl")
cd(dirname(@__FILE__))

# List of open points
"""
- implement linear version of steady-state model (substitution of squared pressures) - Y
- Implement a function to expost calculate linepack for PDE based models - Y

- Harmonize linepack balance (end) condition; especially compare PDE and MISOCP versions - E
- How to incorporate first time step linepack in PDE model - E

- Implement nonconvex with explicit linepack version of model - Y
- update restore_feasibility function - Y
- check SCP as improvement to SLP - Y
- Should we use bus injection vector implementation? - Y
- harmonize rounding?! - E
- harmonize names for Belgian case study folders - E
- Case Study (Large?) - E/Y
- Metrics - E/Y
"""

function store_data(ES; output_file = "training_data_ipopt_QD_linear", append = true, all_pipes = true)
    folder = joinpath(dirname(@__DIR__), "inputs", "neural_network_par", "data", "ipopt")
    #store data
    q_nm = [value(ES.model[:Q_nm][t,p]) for t in ES.T, p in ES.P]
    d = Int(size(q_nm)[1]*size(q_nm)[2])
    q_nm = reshape(q_nm, (d,))
    pr_n = [value(ES.model[:pr][t,ES.pipes[p].start]) for t in ES.T, p in ES.P]
    pr_n = reshape(pr_n, (d,))
    pr_m  = [value(ES.model[:pr][t,ES.pipes[p].stop]) for t in ES.T, p in ES.P]
    pr_m = reshape(pr_m, (d,))
    k_nm = [ES.pipes[p].Knm for t in ES.T, p in ES.P]
    k_nm = reshape(k_nm, (d,))

    t1=hcat(q_nm,pr_n,pr_m,k_nm)
    #concatenate
    #column names: q_nm, pr_n, pr_m, k_nm

    #save to csv
    # if all_pipes
    out_all = joinpath(folder, output_file*"_all_pipes.csv")
    CSV.write(out_all, Tables.table(t1), append = append)
    # else
    tpipes = permutedims(reshape(t1,24,21,4), [1,3,2])
    # tpipes = reshape(t1, :,24,4)
    for p in ES.P
        output_pipe = joinpath(folder, "individual_pipes" , output_file*"_pipe"*string(p)*".csv")
        println("----Storing data for pipe $p----")
        CSV.write(output_pipe, Tables.table(tpipes[:,:,p]), append = append)
    end

end


# gen_data(ES, S = UB, I=2, st=0.01, ths = 0.035, append = true)

########## System must be either 'integrated_power_and_gas' or 'gas_only'
system = "gas_only"

config_dict = Dict(
    ########## Choose case study ##########
    # Gas_tree: 4-node, 3 pipes, 1 source, uniform load #(provided every hour in the input data)
    # Gas_tree_2sources: 4-node, 3 pipes, 2 surce, uniform load #(provided every hour in the input data)
    # Gas_tree_2sources_3min: 4-node, 3 pipes, 2 source, # average load (provided every 3 min in the input data)
    # :case_study => "Gas_tree_compr",
    :case_study => "Belgian", # "IEEE24_Belgium",
    ########## Time discretization ##########
    :timehorizon => 24, # Time horizon (hours)
    :dt => 3600, # Time interval for discretization (seconds)
    ########## Space discretization ##########
    # 0: No space discretization --> discretization segments equal pipeline length
    # 1: Space discretization--> user defines each segments' length
    :space_disc => 0,
    :segment => 10000, #only relevant if space_disc=1
    # Modeltype: choose one of
    # "nonconvex", "MISOCP_bilinear", "MISOCP_McCormick", "SLP", "nonconvex_lp_form", "PWL_incremental",
    :model_type => "nonconvex",
    ########## Choose PDE version ##########
    # Only releveant for model_type nonconvex (1: steady state, 2: quasi-dynamic, 3: transient) 
    :PDE_type => 2,
    ########## Other Parameters ##########
    :print_model => 0, # Control for model printing
    :C_cur_load_gas => 1000, # Value of lost load (kgh)
    :M => 10000, # Parameter for Big-M in MISCOP,
    :obj_type => "linear" #choose between "linear" or "quadratic" objective
)

config_dict_algorithm = Dict{Symbol,Any}()
########## Select additional model parameters for various model formulations
if config_dict[:model_type] == "PWL_incremental"
    ########## Parameters PWL ##########
    config_dict_algorithm[:no_V] = 4 # Number of set points for linearization
elseif config_dict[:model_type] == "SLP"
    ########## Parameters SLP ##########
    config_dict_algorithm[:k_max] = 1000 # Maximum number of iterations
end

# At least the solver name MUST be set
config_dict_solver = Dict(
    :solver_name => "Ipopt",
    :max_cpu_time => 20.0,
)

ES = intialize_energy_system(config_dict, system)
ES.config_dict_solver = config_dict_solver
intialize_data!(ES)

# Random.seed!(1234)
#should probably also store the new demand profile, so that we can compare it to the original one
function gen_data(ES; S=10, i=1, I=10, silent = true, ths = 0.04, UB = 0, st=0.05, output = "training_data_ipopt_QD_linear", append = true)
    S_0 = S
    while true

        
        if S < ths || i > I
            println("#########\nBREAKING\n#########")
            println("S is now $S < $ths, i is now $i great than $I")
            return UB
            break
        end

        intialize_data!(ES)
        println("#########\nITERATION $i:\n#########\n")
        println("----------------\ni = $i, S = $S\n----------------\n\n\n")
        for k in keys(ES.demands)
            d = ES.demands[k]
            # println("$d\n") 
            #sample from uniform distribution, then multiply by 2 in order to get the same average.
            #this changes the load profile randomly for each node, but should, on average, lead to same total demand*S
            d.Qd_gas = d.Qd_gas .* rand(size(d.Qd_gas)[1]) * 2 * S 
        end
        
        run_model!(ES, config_dict_algorithm, silent = silent);
        
        println("Termination status $(termination_status(ES.model))")

        if termination_status(ES.model) != MOI.LOCALLY_SOLVED
            # S=S/2
            println("--------------------------------\nModel not feasible... \nS reduced to $S\n--------------------------------\n")
            UB = S
        end
        if termination_status(ES.model) == MOI.LOCALLY_SOLVED
            println("--------------------------------\nModel is locally feasible! \nStoring data...\n--------------------------------\n")
            #save the data:
            store_data(ES, output_file = output, append = append)
        end
         
        S-=S_0*st #decrease S by 10% of original value for next iteration
        # S+=S*0.1 #increase S by 10% for next iteration
        i+=1
    end

end

# from tests...
# UB = 2.1308126895
# LB = 0.0495
UB = 2.1308126895
UB = 2.5
# gen_data(ES, S = UB, I=450, st=0.01, ths = 0.035)
# println("----------------second run done----------------")
gen_data(ES, S = 2.1, I=450, st=0.01, ths = 0.035)

for i in 1:20
    println("-------------starting run $i----------------")
    gen_data(ES, S = UB, I=450, st=0.01, ths = 0.035)

    println("----------------run $i done----------------")

end
#8040 data points for each pipe
# gen_data(ES, S = UB, I=450, st=0.01, ths = 0.035, output_file)


# UB = 2.1523360499999997
# UB = gen_data(ES, S = UB, I=5, st=0.01)


# UB = gen_data(ES, S = UB, I=5, st=0.01)
# UB
# 2-2*0.1

# function test_iterations(UB, st, T) 
#     for t in 1:T
#          UB-=UB*st
#     end
#     return UB
# end

# test_iterations(UB, 0.01, 350)
