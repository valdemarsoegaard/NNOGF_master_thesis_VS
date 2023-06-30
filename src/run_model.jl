using Pkg
cd(dirname(dirname(@__FILE__)))
Pkg.activate(".") # activate current directory as project directory
# Pkg.add(["JuMP", "Ipopt", "Gurobi", "CSV", "DataFrames", "Statistics", "Plots", "LinearAlgebra", "Dates", "JSON"])
using JuMP, Ipopt, Gurobi, CSV, DataFrames, Statistics, Plots, LinearAlgebra, Dates, JSON, ArgParse, XLSX, Logging
ENV["CPLEX_STUDIO_BINARIES"] = "/apps/cplex/cplex1210/cplex/bin/x86-64_linux"
using CPLEX, HiGHS
include("components.jl")
include("discretization.jl")
include("opt_models_gas.jl")
include("opt_models_integrated_power_and_gas.jl")
include("input_data.jl")
include("algorithms.jl")
include("energy_system.jl")
include("PTDF_matrix.jl")
include("model_functions.jl")
include("output.jl")
cd(dirname(dirname(@__FILE__)))

#decomp files
include("JumpModelToMatrix.jl")
include("DecompMatrix.jl")

using LinearAlgebra, SparseArrays
using ..JumpModelToMatrix, ..DecompMatrix

 


# logger = SimpleLogger(open("log.txt", "w+"))
# global_logger(logger)
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

#ArgParse function to parse command line arguments
s = ArgParseSettings()
@add_arg_table s begin
    "--system"
    help = "System to run\n options: integrated_power_and_gas, gas_only"
    arg_type = String
    default = "gas_only"

    "--case_study"
    help = "Case study to run"
    arg_type = String
    default = "Belgian"# "Belgian"

    "--model_type"
    help = "Model type to run\n options: nonconvex, MISOCP_bilinear, MISOCP_McCormick, SLP, nonconvex_lp_form, PWL_incremental, NN_constrained"
    arg_type = String
    default = "NN_constrained" #nonconvex #NN_constrained

    "--obj_type"
    help = "Objective function type\n options: 'linear': linear objective function, minimizing cost, 2: 'quadratic' objective function, minimizing cost "
    arg_type = String
    default = "linear"

    "--timehorizon"
    help = "Time horizon to run"
    arg_type = Int
    default = 24 #24

    "--dt"
    help = "Time interval for discretization (seconds)"
    arg_type = Int
    default = 3600

    "--space_disc"
    help = "Space discretization\n 0: No space discretization --> discretization segments equal pipeline length\n 1: Space discretization--> user defines each segments' length"
    arg_type = Int
    default = 0

    "--segment"
    help = "only relevant if space_disc=1"
    arg_type = Int
    default = 10000

    "--PDE_type"
    help = "Only releveant for model_type nonconvex (1: steady state, 2: quasi-dynamic, 3: transient)"
    arg_type = Int
    default = 2

    "--print_model"
    help = "Control for model printing"
    arg_type = Int
    default = 0

    "--C_cur_load_gas"
    help = "Value of lost load (kgh)"
    arg_type = Int
    default = 1000

    "--M"
    help = "Parameter for Big-M in MISCOP"
    arg_type = Int
    default = 10000

    "--no_V"
    help = "Number of set points for piecewise linearization"
    arg_type = Int
    default = 4

    "--k_max"
    help = "Maximum number of iterations"
    arg_type = Int
    default = 1000

    "--nn_path"
    help = "Path to neural network model parameters"
    arg_type = String
    default = "best_par_uniform_3_7_v2run_2023_05_05_09_52_58.json" # best_par_linspace_3_7_v2run_2023_05_05_09_28_54#best_par_uniform_3_7_v2run_2023_05_05_09_52_58.json #wandb-summary_run2.json #wandb-summary_ipopt_test.json

    "--solver_name"
    help = "Solver to use\n options: Ipopt, Gurobi"
    arg_type = String
    default = "Gurobi" #Gurobi #Ipopt

    "--nn_M_up"
    help = "upper bound for nn bounds" 
    arg_type = Float64
    default = 100.0

    "--nn_M_low"
    help = "lower bound for nn bounds"
    arg_type = Float64
    default = -100.0

    "--results_path"
    help = "Path to save results"
    arg_type = String
    default = "run_test"

    "--individual_pipes" 
    help = "Use a separate neural network for each pipe"
    arg_type = Bool
    default = false

    "--callback"
    help = "Use callback function"
    arg_type = Bool
    default = false

    "--heuristic_level"
    help = "Level of heuristic to use"
    arg_type = Float64
    default = 0.05

    "--num_threads"
    help = "Number of threads to use"
    arg_type = Int
    default = 1

    "--MIPFocus"
    help = "MIPFocus parameter for Gurobi"
    arg_type = Int
    default = 0

    "--VarBranch"
    help = "VarBranch parameter for Gurobi"
    arg_type = Int
    default = -1

    "--MIPGap"
    help = "MIPGap tolerance for solver"
    arg_type = Float64
    default = 1e-4

    "--presolve"
    help = "presolve parameter for solver"
    arg_type = Int
    default = 1
    
end

parsed_args = parse_args(s,as_symbols = true)

parsed_args[:MIPGap] = 0
# parsed_args[:callback] = true
# parsed_args[:num_threads] = 6
# parsed_args[:heuristic_level] = 0.8
# parsed_args[:nn_path] = "best_par_linspace_3_7_all_v3run_2023_04_21_16_56_44.json"
# parsed_args[:nn_path] = "best_par_uniform_toy_examplerun_2023_05_26_11_17_50.json"
parsed_args[:model_type] = "nonconvex" #" nonconvex, MISOCP_bilinear, MISOCP_McCormick, SLP, nonconvex_lp_form, PWL_incremental, NN_constrained"
# parsed_args[:solver_name] = "Gurobi" #"Ipopt, Gurobi, HiGHS
# parsed_args[:solver_name] = "CPLEX"
parsed_args[:solver_name] = "Ipopt"

parsed_args[:case_study] = "Gas_tree_compr" #"gas_tree_compr"
# parsed_args[:results_path] = "$(parsed_args[:solver_name][1])_case_$(parsed_args[:case_study])_mt_$(parsed_args[:model_type])_BP$(parsed_args[:no_V])"
# parsed_args[:results_path] = "$(parsed_args[:solver_name][1])_case_$(parsed_args[:case_study])_mt_$(parsed_args[:model_type])"
# parsed_args[:results_path] = "$(parsed_args[:solver_name][1])_case_$(parsed_args[:case_study])_CB_$(parsed_args[:callback] == true ? "t" : "f")_mt_$(parsed_args[:model_type])_ps$(parsed_args[:presolve])"
# parsed_args[:results_path] = "test123"
# parsed_args[:dt] = 900
# parsed_args[:results_path]  = "run_test_ipopt_nonconvex"
# parsed_args[:nn_path] = "best_par_combined_3_7run_2023_05_04_16_37_16.json"

# println("Parsed arguments: ")
# print(parsed_args)

if parsed_args[:callback] == true
    set_logger(parsed_args[:results_path])
end

@info "Starting run"
@info "Parsed arguments: "
@info parsed_args

########## System must be either 'integrated_power_and_gas' or 'gas_only'
# system = "gas_only"

config_dict = Dict(
    ########## Choose case study ##########
    # Gas_tree: 4-node, 3 pipes, 1 source, uniform load #(provided every hour in the input data)
    # Gas_tree_2sources: 4-node, 3 pipes, 2 surce, uniform load #(provided every hour in the input data)
    # Gas_tree_2sources_3min: 4-node, 3 pipes, 2 source, # average load (provided every 3 min in the input data)
    # :case_study => "Gas_tree_compr",
    :case_study => parsed_args[:case_study], # "IEEE24_Belgium",
    ########## Time discretization ##########
    # :timehorizon => 24, # Time horizon (hours)
    :timehorizon => parsed_args[:timehorizon], # Time horizon (hours)
    :dt => parsed_args[:dt], # Time interval for discretization (seconds)
    ########## Space discretization ##########
    # 0: No space discretization --> discretization segments equal pipeline length
    # 1: Space discretization--> user defines each segments' length
    :space_disc => parsed_args[:space_disc],
    :segment => parsed_args[:segment], #only relevant if space_disc=1
    # Modeltype: choose one of
    # "nonconvex", "MISOCP_bilinear", "MISOCP_McCormick", "SLP", "nonconvex_lp_form", "PWL_incremental", "NN_constrained"
    :model_type => parsed_args[:model_type],
    # :model_type => "PWL_incremental",
    ########## Choose PDE version ##########
    # Only releveant for model_type nonconvex (1: steady state, 2: quasi-dynamic, 3: transient) 
    :PDE_type => parsed_args[:PDE_type],
    ########## Other Parameters ##########
    :print_model => parsed_args[:print_model], # Control for model printing
    :C_cur_load_gas => parsed_args[:C_cur_load_gas], # Value of lost load (kgh)
    :M => parsed_args[:M], # Parameter for Big-M in MISCOP
    :obj_type => parsed_args[:obj_type], #choose between "linear" or "quadratic" objective
)

config_dict_algorithm = Dict{Symbol,Any}()
config_dict_algorithm[:callback] = parsed_args[:callback]
config_dict_algorithm[:heuristic_level] = parsed_args[:heuristic_level]
config_dict_algorithm[:nn_M_up] = parsed_args[:nn_M_up]
config_dict_algorithm[:nn_M_low] = parsed_args[:nn_M_low]
config_dict_algorithm[:individual_pipes] = parsed_args[:individual_pipes]
########## Select additional model parameters for various model formulations
if config_dict[:model_type] == "PWL_incremental"
    ########## Parameters PWL ##########
    config_dict_algorithm[:no_V] = parsed_args[:no_V] # Number of set points for linearization
elseif config_dict[:model_type] == "SLP"
    ########## Parameters SLP ##########
    config_dict_algorithm[:k_max] = 1000 # Maximum number of iterations
elseif config_dict[:model_type] == "NN_constrained"
    ########## Parameters SLP ##########
    nn_path = parsed_args[:nn_path] #change this to the path of the wandb-summary file
    config_dict_algorithm[:NN_path] = joinpath(dirname(@__DIR__), "inputs", "neural_network_par","parameters", nn_path)
end

# At least the solver name MUST be set
config_dict_solver = Dict{Symbol, Any}()
config_dict_solver[:solver_name] = parsed_args[:solver_name]
# config_dict_solver = Dict(
#     :solver_name => parsed_args[:solver_name],#Ipopt
#     # :max_cpu_time => 100.0, #only run with Ipopt...
#     # :NonConvex => 2,
#     # :TimeLimit => 100,
#     # :FeasibilityTol => 1e-9
# )
if parsed_args[:solver_name] == "Gurobi"
    config_dict_solver[:MIPFocus] = parsed_args[:MIPFocus]
    config_dict_solver[:VarBranch] = parsed_args[:VarBranch]
    config_dict_solver[:Heuristics] = parsed_args[:heuristic_level]
    config_dict_solver[:Threads] = parsed_args[:num_threads]
    config_dict_solver[:MIPGap] = parsed_args[:MIPGap]

    if parsed_args[:presolve] == 0
        config_dict_solver[:Presolve] = parsed_args[:presolve]
    else
        config_dict_solver[:Presolve] = -1 # if "on" use the dfault value
    end
end
if parsed_args[:solver_name] == "CPLEX"
    config_dict_solver[:CPXPARAM_Threads] = parsed_args[:num_threads]
    config_dict_solver[:CPXPARAM_MIP_Tolerances_MIPGap] = parsed_args[:MIPGap]
    config_dict_solver[:CPXPARAM_Preprocessing_Presolve] = parsed_args[:presolve]

    config_dict_solver[:CPXPARAM_MIP_Display] = 5
    config_dict_solver[:CPXPARAM_MIP_Strategy_KappaStats] = 0
end

ES = intialize_energy_system(config_dict, parsed_args[:system])
ES.config_dict_solver = config_dict_solver

intialize_data!(ES)
# # config_dict_algorithm[:individual_pipes] = true
# m = build_NN_constrained_model(ES, config_dict_algorithm[:NN_path], 
#                                                 M_up = config_dict_algorithm[:nn_M_up], 
#                                                 M_low = config_dict_algorithm[:nn_M_low],
#                                                 individual_pipes = config_dict_algorithm[:individual_pipes])
# run_model!(ES, config_dict_algorithm, LP_relax = true); 
run_model!(ES, config_dict_algorithm);
# run_model!(ES, config_dict_algorithm, var_fix = true);

for i in 1:21
    println(ES.pipes[i].S)
end




# index(ES.model[:b][bid])
# index(ES.model[:b])
# bid= CartesianIndex(11, 20, 2, 5)
# bid= CartesianIndex(1, 3, 2, 5)
# index(ES.model[:b][bid])

# xn_val = value.(ES.model[:xn])
# xp_val = value.(ES.model[:xp])
# y_val = value.(ES.model[:y])
# y_hat_val = value.(ES.model[:y_hat])
# b_val = value.(ES.model[:b])

# ReLU_mask_pos = xn_val .== 0 .&& xp_val .> 0; # if the negative part is zero, then the ReLU is not violated, as the positive part must be correct. Maybe do between +/- 1e-6 instead of zero?
# sum(ReLU_mask_pos)
# #Find ReLUs that are not violated where positive part is zero and negative part is not zero
# ReLU_mask_neg = xp_val .== 0 .&& xn_val .< 0;
# sum(ReLU_mask_neg)
# # find ReLUs that are not violated and both positive and negative parts are zero
# ReLU_mask_both = xp_val .== xn_val .== 0;
# sum(ReLU_mask_both)

# sum(ReLU_mask_pos) + sum(ReLU_mask_neg) + sum(ReLU_mask_both)
# #filter fractional binaries corresponding to these ReLUs
# b_mask_pos = 0 .< b_val[ReLU_mask_pos] .< 1;
# sum(b_mask_pos)
# # println(sum(b_mask_pos))
# # @info "b_mask_pos: $(sum(b_mask_pos))"
# b_mask_neg = 0 .< b_val[ReLU_mask_neg] .< 1;
# sum(b_mask_neg)
# # println(sum(b_mask_neg))
# # @info "b_mask_neg: $(sum(b_mask_neg))"
# b_mask_both = 0 .< b_val[ReLU_mask_both] .< 1;
# sum(b_mask_both)
# sum(b_mask_pos) + sum(b_mask_neg) + sum(b_mask_both)

# #get the indices of the fractional binaries
# #get the variables corresponding to these indices
# b_frac_pos = m[:b][ReLU_mask_pos][b_mask_pos]
# b_frac_neg = m[:b][ReLU_mask_neg][b_mask_neg]
# b_frac_both = m[:b][ReLU_mask_both][b_mask_both]
# #MAKE SURE THAT INDEX DOES NOT APPEAR TWICE, WHEN WE HAVE BOTH POSITIVE AND NEGATIVE RELU PARTS = 0 
# @info "Fractional ReLUs: positive $(sum(b_mask_pos)), negative $(sum(b_mask_neg)) and both $(sum(b_mask_both)) "
# # b_var_sel = [b_frac_pos; b_frac_neg] #
# b_var_sel = [b_frac_pos; b_frac_neg; b_frac_both]
# # b_var_val = [floor.(Int,b_val[ReLU_mask_pos][b_mask_pos]); ceil.(Int,b_val[ReLU_mask_neg][b_mask_neg])] #might need to be float instead of int
# b_var_val = [b_val[ReLU_mask_pos][b_mask_pos]; b_val[ReLU_mask_neg][b_mask_neg];b_val[ReLU_mask_both][b_mask_both]]

# sum(0.0 .< b_val .< 1.0)

# minimum(b_val[0.0 .< b_val .< 1.0])
# maximum(b_val[0.0 .< b_val .< 1.0])

# L = 0.99991
# L = 0.1
# b_val[0.0 .< b_val .< 1.0][L .> b_val[0.0 .< b_val .< 1.0]]
# y_hat_val[0.0 .< b_val .< 1.0][L .> b_val[0.0 .< b_val .< 1.0]]
# y_val[0.0 .< b_val .< 1.0][L .> b_val[0.0 .< b_val .< 1.0]]




# minimum(y_hat_val[0.0 .< b_val .< 1.0])
# minimum(y_val[0.0 .< b_val .< 1.0])
# minimum(xp_val[0.0 .< b_val .< 1.0])
# minimum(xn_val[0.0 .< b_val .< 1.0])


# argmax(b_val[0.0 .< b_val .< 1.0])

# b_val[0.0 .< b_val .< 1.0][2306]
# y_val[0.0 .< b_val .< 1.0][2306]
# y_hat_val[0.0 .< b_val .< 1.0][2306]

# sum(1.0 .< d_fix)
# sum(d_fix .< 1.0)


# b_frac = 0.0 .< b_val .< 1.0 


# # y_hat_val_p = min.(y_hat_val, 0)
# # y_hat_val_p = min.(y_hat_val, 0)
# # d = (.-xn_val.+xp_val.+1)./(y_val.+1)
# y_hat_val_p = min.(y_hat_val, 0)
# # d_fix = (.-xn_val.+xp_val.+1)./(y_val-y_hat_val_p.+1)
# d_fix = (.-xn_val[b_frac].+xp_val[b_frac].+1)./(y_val[b_frac]-y_hat_val_p[b_frac].+1)
# i = argmax(d_fix)
# i_ = findall(d_fix .== 1.0)

# m[:b][b_frac][i]

# d_fix[i_]
# b_val[i_]

# d_fix
# # d[:,:,2,5]
# d_fix[:,:,2,5]

# y_val[1,7,end,end]
# y_hat_val[1,7,end,end]#[end,end,end,end]
# xn_val[1,7,end,end]#[end,end,end,end]
# xp_val[1,7,end,end]#[end,end,end,end]
# d_fix[1,7,end,end]



# d_id = argmax(d_fix)
# d[d_id]
# y_val[d_id]
# xn_val[d_id]
# xp_val[d_id]
# y_hat_val[d_id]

#----------------- PLOT RESULTS -----------------#
println(parsed_args[:results_path])
# parsed_args[:results_path] = "uniform_data_CPLEX_case_Gas_tree_compr_CB_false_modeltype_NN_constrained_var_fix"
Plotres(ES.N,ES.pr,ES.N_dt, parsed_args[:results_path], save = true)
save_results(ES, parsed_args[:results_path])

#----------------- RESTORE FEASIBILITY -----------------#
# attribute_dict = Dict(
#             :solver_name => "Ipopt",
#             :max_cpu_time => 600.0,
#             :max_iter => 10000,
#         )

# ES_res = intialize_energy_system(config_dict, parsed_args[:system])
# ES_res.config_dict_solver = config_dict_solver
# intialize_data!(ES_res)
# ES_res = restore_feasibility(ES_res, ES.Q_w, ES.Q_nm, ES.pr, attribute_dict = attribute_dict);

# # ES_res = restore_feasibility(ES, attribute_dict = attribute_dict);

# if termination_status(ES_res.model) != MOI.LOCALLY_SOLVED
#     # println("No feasible solution found")
#     @info "No feasible solution found"
# else
#     Plotres(ES_res.N,ES_res.pr,ES_res.N_dt, parsed_args[:results_path], save = true, name = "res_feas_nodal_pres")
#     save_results(ES_res, parsed_args[:results_path], "res_feas_results.xlsx")
# end

#------------------ OTHER STUFF------------------------------#

# a = 2

#ADD CONVERGENCE PLOT!!
# Building constraint Matrix

# m = build_NN_constrained_model(ES)
# # A, b, vecObj, vecVarLB, vecVarUB, vecIsInt, sense, varNames = getConstraintMatrix(m)
# # CSV.write("Am.csv", Tables.table(Matrix(A)))
# nnogf, nndic = getConstraintMatrix(m)
# A  = Matrix(nnogf.A)
# CSV.write("Am_small_t1_named.csv", Tables.table(A, header = nnogf.varNames))

# #Calculate linepack (ex-post) for quasi-dynamic and dynamic model with PDE formulation
# if (var>1)&&(lp_form==0)
#     LP=[round(pipes[p].LL/pipes[p].V_m*value(pr_avg[t,p]), digits=2) for t in T, p in P] #kg
# end

## Create dictionary with all results

## Save results for post processing
# using FileIO, JLD2
#name_file=string("results_", case_study,"_var",var,"_averageload","_dt",dt,".jld2")

#FileIO.save(name_file, "case_study",case_study,"dt",dt,"space_disc",space_disc,"lp_form",lp_form,"obj_val", obj_val, "solve_time",solve_time,"Q_w", Q_w,"pr", pr, "pr_avg",pr_avg,"Q_nm", Q_nm, "Q_nm_in",Q_nm_in,"Q_nm_out", Q_nm_out, "Q_c", Q_c, "Qd_cur_gas", Qd_cur_gas, "LP",LP )
#FIleIO.save(name_file,"dt",dt,"space_disc",space_disc,case_study",case_study,"lp_form",lp_form,"obj_val",obj_val, "solve_time",solve_time, "Q_w", Q_w,"pr", pr,"Q_nm", Q_nm, "Q_nm_in",Q_nm_in,"Q_nm_out", Q_nm_out, "pr_avg",pr_avg,"Qd_cur_gas", Qd_cur_gas, "Q_c", Q_c,"LP",LP )
#a = 1
#FileIO.save(name_file,"a",a)
#b = FileIO.load(name_file,"a") 

#-----------------Testing NN for each pipe

# "C:\Users\valde\OneDrive - Danmarks Tekniske Universitet\GitHub\OGF\inputs\neural_network_par\parameters\linspace\best_par_linspace_3_7_v2run_2023_05_05_09_28_54.json"

# config_dict_algorithm[:NN_path] = joinpath(dirname(@__DIR__), "inputs", "neural_network_par","parameters", "best_par_linspace_3_7_v2run_2023_05_05_09_28_54.json")
# config_dict_algorithm[:individual_pipes] = false
# build_NN_constrained_model(ES, config_dict_algorithm[:NN_path], 
#                                                 M_up = config_dict_algorithm[:nn_M_up], 
#                                                 M_low = config_dict_algorithm[:nn_M_low],
#                                                 individual_pipes = config_dict_algorithm[:individual_pipes])
                                               

# config_dict_algorithm[:individual_pipes] = true
# build_NN_constrained_model(ES, config_dict_algorithm[:NN_path], 
#                                                 M_up = config_dict_algorithm[:nn_M_up], 
#                                                 M_low = config_dict_algorithm[:nn_M_low],
#                                                 individual_pipes = config_dict_algorithm[:individual_pipes])


#----------------- Plot a subset of pipes---------------------------------
# N = ES,pr,N_dt
# folder = joinpath(dirname(@__DIR__), "output", output_folder)

# fig_name =  "run_test_gurobi_nonconvex"
# #clrs = [:red :blue :green  :cyan :magenta]
# fig = plot(title="Nodal pressure", xlabel="Time step", ylabel="Pressure [MPa]", legend=:outerright, dpi = 300)
#     for n in (6,7,4)
#         #plot!(1:N_dt, pr[:,n], lc=clrs[n], label="Node $n")
#         plot!(1:ES.N_dt, ES.pr[:,n], label="Node $n")
#     end
#     display(fig)
#     if save
#         savefig(fig, fig_name) #how to set dpi higher
#     end
#     function Plotres(N,pr,N_dt, output_folder :: String; save = false, name = "Nodal_pressure.png") #add more plots?

#----------------- bounds....---------------------------------
# m = build_NN_constrained_model(ES)
#IPOPT_test
# config_dict_algorithm[:nn_M_up] = 1.95
# config_dict_algorithm[:nn_M_low] = -1.73
# uniform...
# config_dict_algorithm[:nn_M_up] = 0.92025
# config_dict_algorithm[:nn_M_low] = -0.99995


#------------- Branching / Branching Rules ----------------

# b = value.(ES.model[:b])
# y = value.(ES.model[:y])
# y_hat = value.(ES.model[:y_hat])
# xn = value.(ES.model[:xn])
# xp = value.(ES.model[:xp])
# y0 = value.(ES.model[:y0])

# b[:,:,1,1]
# y_hat[:,:,1,1]
# y[:,:,1,1]
# xn[:,:,1,1]
# xp[:,:,1,1]

# xn[:,:,1,1]+xp[:,:,1,1]

# y_hat[:,:,1,1].>0
# y_hat[:,:,1,1]

# yp = max.(0,y_hat[:,:,1,1])
# 15,8

# y_hat[:,:,1,1]
# y[:,:,1,1]
# y[:,:,1,1]-yp

# y_hat[15,8,1,1]
# xn[15,8,1,1]
# xp[15,8,1,1]
# y[15,8,1,1]
# b[15,8,1,1]

# y_hat[23,1,1,1]
# xn[23,1,1,1]
# xp[23,1,1,1]
# y[23,1,1,1]
# b[23,1,1,1]

# b
# round.(b,6) .== 0.517279

# b[20,10,2,5]
# y_hat[20,10,2,5]
# y[20,10,2,5]
# xn[20,10,2,5]
# xp[20,10,2,5]


# println("Example of bad branching:")
# println("b[20,10,2,5] = ", b[20,10,2,5])
# println("y_hat[20,10,2,5] = ", y_hat[20,10,2,5])
# println("y[20,10,2,5] = ", y[20,10,2,5])
# println("xn[20,10,2,5] = ", xn[20,10,2,5])
# println("xp[20,10,2,5] = ", xp[20,10,2,5])
# println("Branching rule (xp - xn)/y : $((-xn[20,10,2,5]+xp[20,10,2,5])/y[20,10,2,5])")
# println("Branching rule (xp - xn)-y : $((-xn[20,10,2,5]+xp[20,10,2,5])-y[20,10,2,5])")
# println("Branching rule (xp - xn) : $((-xn[20,10,2,5]+xp[20,10,2,5]))")

# #yh =-0.5, y = 0.2, xn = -0.7,xp = 0.2
# println("This would also be valid for b = 1, thus it makes no sense to branch on this variable first (this mapping is already feasible)")

# println("Example of better branching:")
# println("b[23,1,1,1] = ", b[23,1,1,1])
# println("y_hat[23,1,1,1] = ", y_hat[23,1,1,1])
# println("y[23,1,1,1] = ", y[23,1,1,1])
# println("xn[23,1,1,1] = ", xn[23,1,1,1])
# println("xp[23,1,1,1] = ", xp[23,1,1,1])
# println("Branching rule (xp - xn)/y : $((-xn[23,1,1,1]+xp[23,1,1,1])/y[23,1,1,1])")
# println("Branching rule (xp - xn)-y : $((-xn[23,1,1,1]+xp[23,1,1,1])-y[23,1,1,1])")
# println("Branching rule (xp - xn) : $((-xn[23,1,1,1]+xp[23,1,1,1]))")
# println("This is not valid for b = 1, thus it makes sense to branch on this variable first")


# # y_hat[23,1,1,1]
# # xn[23,1,1,1]
# # xp[23,1,1,1]
# # y[23,1,1,1]
# # b[23,1,1,1]

# println("Another example of good branching:")
# println("b[15,8,1,1] = ", b[15,8,1,1])
# println("y_hat[15,8,1,1] = ", y_hat[15,8,1,1])
# println("y[15,8,1,1] = ", y[15,8,1,1])
# println("xn[15,8,1,1] = ", xn[15,8,1,1])
# println("xp[15,8,1,1] = ", xp[15,8,1,1])
# println("Branching rule (xp - xn)/y : $((-xn[15,8,1,1]+xp[15,8,1,1])/y[15,8,1,1])")
# println("Branching rule (xp - xn)-y : $((-xn[15,8,1,1]+xp[15,8,1,1])-y[15,8,1,1])")
# println("Branching rule (xp - xn) : $((-xn[15,8,1,1]+xp[15,8,1,1]))")
# println("This is not valid for b = 1, thus it makes sense to branch on this variable first")


# scatter([y_hat[15,8,1,1]],[y[15,8,1,1]])
# display()
# y_hat[15,8,1,1]
# xn[15,8,1,1]
# xp[15,8,1,1]
# y[15,8,1,1]
# b[15,8,1,1]

# (-xn[15,8,1,1]+xp[15,8,1,1])/y[15,8,1,1]
# (-xn[23,1,1,1]+xp[23,1,1,1])/y[23,1,1,1]

# #---------- Counting Violated ReLUs-------------------
# yp_all = max.(0,y_hat);

# y_wrong = y - yp_all

# count_rel_wrong = y_wrong.>0
# count_bin_wrong = 1 .> b .> 0


# println(" Number of wrong ReLUs: $(sum(count_rel_wrong))")
# b_val = value.(ES.model[:b])
# y_val = value.(ES.model[:y])
# y_hat_val = value.(ES.model[:y_hat])
# xn_val = value.(ES.model[:xn])
# xp_val = value.(ES.model[:xp])
# y0_val = value.(ES.model[:y0])


# b = ES.model[:b]
# y = ES.model[:y]
# y_hat = ES.model[:y_hat]
# xn = ES.model[:xn]
# xp = ES.model[:xp]
# y0 = ES.model[:y0]

# prod(size(y_hat .- xp))


# ReLU_mask_pos = y_val - xp_val - xn_val .== 0
# sum(ReLU_mask_pos)
# y[ReLU_mask_pos]
# y_hat[ReLU_mask_pos]
#compute mask for violated b with positive y_hat
# b_mask_pos = 0 .< b_val[ReLU_mask_pos] .< 1
# sum(b_mask_pos)
# b[ReLU_mask_pos][b_mask_pos]

# b_frac_pos = b[ReLU_mask_pos][b_mask_pos] 

# ReLU_mask_neg = xp_val .== 0
# # sum(ReLU_mask_neg)
# # y[ReLU_mask_neg]
# # y_hat[ReLU_mask_neg]
# b_mask_neg = 0 .< b_val[ReLU_mask_neg] .< 1
# # sum(b_mask_neg)
# # b[ReLU_mask_neg][b_mask_neg]

# b_frac_neg = b_val[ReLU_mask_neg][b_mask_neg]

# if isempty(b_frac_neg)
#     b_var_sel = b_frac_pos
#     b_var_val = [floor.(Int,b_val[ReLU_mask_pos][b_mask_pos])]
# else
#     b_var_sel = [b_frac_pos b_frac_neg]
#     b_var_val = [floor.(Int,b_val[ReLU_mask_pos][b_mask_pos]) ceil.(Int,b_val[ReLU_mask_neg][b_mask_neg])]
# end

# [b_frac_pos; b_frac_neg]
# [floor.(Int,b[ReLU_mask_pos][b_mask_pos]); ceil.(Int,b[ReLU_mask_pos][b_mask_pos])]


# ReLU_mask = y - xp .== 0
# sum(ReLU_mask)



# mask
# b_mask = 0 .< b[ReLU_mask] .< 1

# xp[ReLU_mask]
# xn[ReLU_mask]
# y[ReLU_mask]
# y_hat[ReLU_mask]


# b[ReLU_mask]

# b[ReLU_mask][b_mask]

# y[ReLU_mask][b_mask]

# p_mask = xp[ReLU_mask][b_mask] .> 0

# n_mask = xn[ReLU_mask][b_mask] .< 0

# sum(n_mask)
# sum(b[ReLU_mask][b_mask][p_mask])


# y[:,end,end,end]
# y_hat[:,end,end,end]
# sum(count_bin_wrong)

# prod(size(b))
# prod(size(count_rel_wrong))

# for k in keys(ES.pipes)
#     println("Pipe $k, Knm value:$(ES.pipes[k].Knm)")
# end


# for k in 1:21
#     println("Pipe $k, Knm value:$(ES.pipes[k].Knm)")
# end

# for k in 1:21
#     println(ES.pipes[k].Knm)
# end