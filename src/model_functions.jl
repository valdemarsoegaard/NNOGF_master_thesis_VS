include("energy_system.jl")

function run_model!(ES::EnergySystem, config_dict_algorithm; model_type=nothing, silent = false, LP_relax = false, var_fix = false)
    if isnothing(model_type)
        model_type = ES.model_type
    end

    if model_type in ["SLP"] # Add other algorithms if exist
        ES = SLP(ES, config_dict_algorithm[:k_max])
    else
        if model_type == "nonconvex"
            m = build_PDE_model(ES)
        # elseif model_type == "nonconvex_lp_form"
        #     ES.obj_val, ES.solve_time, ES.Q_w, ES.Q_nm, ES.Q_nm_in, ES.Q_nm_out,
        #         ES.Qd_cur_gas, ES.Q_c, ES.pr, ES.LP = nonconvex_linepack(ES)
        elseif model_type == "MISOCP_bilinear"
            m = build_MISOCP_model(ES)
        elseif model_type == "PWL_incremental"
            m = build_PWL_incremental_model(ES, config_dict_algorithm[:no_V])
        elseif model_type == "MISOCP_McCormick"
            m = build_MISOCP_McCormick_model(ES)
        elseif model_type == "NN_constrained"
            m = build_NN_constrained_model(ES, config_dict_algorithm[:NN_path], 
                                                M_up = config_dict_algorithm[:nn_M_up], 
                                                M_low = config_dict_algorithm[:nn_M_low],
                                                individual_pipes = config_dict_algorithm[:individual_pipes],
                                                heuristic_level = config_dict_algorithm[:heuristic_level])
        else 
            throw(ValueError("Passed model_type is $(ES.model_type). Must be one of
            'nonconvex', 'MISOCP_bilinear', 'MISOCP_McCormick', 'nonconvex_lp_form', 'PWL_incremental', 'NN_constrained'."  ))
        end
        if silent
            MOI.set(m, MOI.Silent(), true)
        end
        undo = nothing
        if LP_relax
            undo = relax_integrality(m)
        end
        if config_dict_algorithm[:callback]
            # println("Callback function is used")
            @info "Callback function is used"
            if solver_name(m) == "Gurobi" #check the syntax... 
                # println("Gurobi optimizer")
                @info "Gurobi optimizer"
                MOI.set(m, Gurobi.CallbackFunction(), my_callback_function)

            elseif solver_name(m) == "CPLEX"
                # println("CPLEX optimizer")
                @info "CPLEX optimizer"
                # set_optimizer(m, CPLEX.Optimizer)
                # cb_calls = Clong[]
                
                # MOI.set(m, MOI.NumberOfThreads(), 1)
                set_optimizer_attribute(m, "CPXPARAM_Threads", 1) #need to run single thread when using callbacks
                if var_fix
                    @info "callback w/ variable fixing is used"
                    MOI.set(m, CPLEX.CallbackFunction(), my_callback_function_vf)
                else
                    @info "Callback function is used"
                    MOI.set(m, CPLEX.CallbackFunction(), my_callback_function)
                end
                
            else
                throw(ArgumentError("Callback function is not implemented for this solver"))
            end
        end
        global m #This is needed for the callback... 
        
        solve_model!(m)
        ES = write_solution(ES, m, model_type);
    end

    return ES, undo
end


function solve_model!(m)
    @time optimize!(m)
    status = termination_status(m)
    # println(status)
    # println(raw_status(m))
    @info status
    @info raw_status(m)
    return
end

function write_solution_gas(ES::EnergySystem, m, model_type)
    # Model statistics
    ES.obj_val = objective_value(m)
    ES.solve_time =  JuMP.solve_time(m)

    # Common variables in all models
    ES.Q_w = [value(m[:Q_w][t,s]) for t in ES.T, s in ES.S]
    ES.pr = [value(m[:pr][t,n]) for t in ES.T, n in ES.N]
    ES.Qd_cur_gas = [value(m[:Qd_cur_gas][t,d]) for t in ES.T, d in ES.D]
    ES.Q_c = isempty(ES.compressors)==true ? [] : [
        value(m[:Q_c][t,c]) for t in ES.T, c in ES.C]
    ES.model = m
    
    # Model-specific variables
    if model_type in ["MISOCP_bilinear", "MISOCP_McCormick"]

        ES.LP = [value(m[:LP][t,p]) for t in ES.T, p in ES.P]
        ES.y_nm = [value(m[:y_nm][t,p]) for t in ES.T, p in ES.P]
        ES.ϕ⁺ = [value(m[:ϕ⁺][t,p]) for t in ES.T, p in ES.P]
        ES.ϕ⁻ = [value(m[:ϕ⁻][t,p]) for t in ES.T, p in ES.P]
        ES.Q_nm = [value(m[:Q_nm][t,p]) for t in ES.T, p in ES.P]
        ES.Q_nm⁺ = [value(m[:Q_nm⁺][t,p]) for t in ES.T, p in ES.P]
        ES.Q_nm_in⁺ = [(value(m[:Q_nm_in⁺][t,p])) for t in ES.T, p in ES.P]
        ES.Q_nm_out⁺ = [(value(m[:Q_nm_out⁺][t,p])) for t in ES.T, p in ES.P]
        ES.Q_nm⁻ = [(value(m[:Q_nm⁻][t,p])) for t in ES.T, p in ES.P]
        ES.Q_nm_in⁻ = [(value(m[:Q_nm_in⁻][t,p])) for t in ES.T, p in ES.P]
        ES.Q_nm_out⁻ = [(value(m[:Q_nm_out⁻][t,p])) for t in ES.T, p in ES.P]
        out_dic = [:Q_w, :pr, :Qd_cur_gas,  :LP, :y_nm, :ϕ⁺, :ϕ⁻, :Q_nm, :Q_nm⁺, :Q_nm_in⁺,
         :Q_nm_out⁺, :Q_nm⁻, :Q_nm_in⁻, :Q_nm_out⁻] 
    elseif model_type in ["PWL_incremental"]

        ES.Q_nm = [(value(m[:Q_nm][t,p])) for t in ES.T, p in ES.P]
        ES.Q_nm_in = [(value(m[:Q_nm_in][t,p])) for t in ES.T, p in ES.P]
        ES.Q_nm_out = [(value(m[:Q_nm_out][t,p])) for t in ES.T, p in ES.P]
        ES.LP = [value(m[:LP][t,p]) for t in ES.T, p in ES.P]
        out_dic = [:Q_w, :pr, :Qd_cur_gas,  :Q_nm, :Q_nm_in, :Q_nm_out, :LP]
    elseif model_type in ["nonconvex", "SLP"]

        ES.Q_nm = [(value(m[:Q_nm][t,p])) for t in ES.T, p in ES.P]
        ES.Q_nm_in = [(value(m[:Q_nm_in][t,p])) for t in ES.T, p in ES.P]
        ES.Q_nm_out = [(value(m[:Q_nm_out][t,p])) for t in ES.T, p in ES.P]
        ES.pr_avg = [(value(m[:pr_avg][t,p])) for t in ES.T, p in ES.P]
        out_dic = [:Q_w, :pr, :Qd_cur_gas,  :Q_nm, :Q_nm_in, :Q_nm_out, :pr_avg]
        
    elseif model_type in ["NN_constrained"]

        ES.Q_nm = [(value(m[:Q_nm][t,p])) for t in ES.T, p in ES.P]
        ES.Q_nm_in = [(value(m[:Q_nm_in][t,p])) for t in ES.T, p in ES.P]
        ES.Q_nm_out = [(value(m[:Q_nm_out][t,p])) for t in ES.T, p in ES.P]
        ES.LP = [value(m[:LP][t,p]) for t in ES.T, p in ES.P]
        out_dic = [:Q_w, :pr, :Qd_cur_gas,  :Q_nm, :Q_nm_in, :Q_nm_out, :LP]#maybe add binary values

    end
    
    if !isempty(ES.compressors)
        push!(out_dic, :Q_c)
    end 

    ES.out_dic = out_dic
    return ES

end


function write_solution_power(ES::IntegratedEnergySystem, m, model_type)
    """
    Extracts and writes model solutions for the power system side on the
    IntegratedEnergySystem object.
    """

    GFPPs = [g for g in ES.G if ES.generators[g].type == "NGFPP"]
    ES.Qd_gfpp = [(value(ES.generators[gg].eta*m[:p][gg, t])) for t in ES.T, gg in GFPPs]
    ES.p = [(value(m[:p][g,t])) for t in ES.T, g in ES.G]
    ES.w = [(value(m[:w][j,t])) for t in ES.T, j in ES.J]

    return ES
end


function write_solution(ES::EnergySystem, m, model_type)
    ES = write_solution_gas(ES, m, model_type)
    if typeof(ES) == IntegratedEnergySystem
        ES = write_solution_power(ES, m, model_type)
    end
    return ES
end

function set_solver(m::JuMP.Model, ES::EnergySystem; attribute_dict::Dict=Dict())
    """
    Sets solver for the optimization model.
    """

    if isempty(attribute_dict)
        attribute_dict = copy(ES.config_dict_solver)
    end
    solver_name = attribute_dict[:solver_name]

    if solver_name == "Ipopt"
        set_optimizer(m, Ipopt.Optimizer)
    elseif solver_name == "Gurobi"
        m = direct_model(Gurobi.Optimizer())
        if ES.model_type in ["nonconvex", "MISOCP_bilinear"]
            set_optimizer_attribute(m, "NonConvex", 2)
        end
        # set_optimizer_attribute(m, "Threads", attribute_dict[:num_threads])
    elseif solver_name == "HiGHS"
        # set_optimizer(m, HiGHS.Optimizer)
        m = direct_model(HiGHS.Optimizer())
        elseif solver_name == "CPLEX"
        # set_optimizer(m, CPLEX.Optimizer)
        CPLX_ENV = CPLEX.Env()
        m = direct_model(CPLEX.Optimizer(CPLX_ENV))
        # status = CPXsetlogfilename(CPLX_ENV, "logfile.txt", "w");
        # set_optimizer_attribute(m, "CPXPARAM_Emphasis_MIP", 1) #Emphasise feasible solutions
        # set_optimizer_attribute(m, "Heuristics", heuristic_level)
    else
        throw(ArgumentError("Solver must be either 'Ipopt', 'Gurobi', 'HiGHS' or 'CPLEX'."))
    end
    delete!(attribute_dict, :solver_name)
    for (key, val) in attribute_dict
        print(key, " ", val, "\n")
        set_optimizer_attribute(m, String(key), val)
    end

    return m 
end


function restore_feasibility(ES::EnergySystem, Qw, Qnm, pr; attribute_dict::Dict=Dict())
    
    m = build_PDE_model(ES)
    for t in ES.T
        for n in ES.N set_start_value(m[:pr][t,n], pr[t,n]) end
        # for s in ES.S set_start_value(m[:Q_w][t,s], Qw[t,s]) end
        # for p in ES.P set_start_value(m[:Q_nm][t,p], Qnm[t,p]) end #try only pressure... 
    end

    if isempty(attribute_dict)
        attribute_dict = Dict(
            :solver_name => "Ipopt",
            :max_cpu_time => 600.0,
        )
    end
    set_solver(m, ES; attribute_dict=attribute_dict)

    solve_model!(m)
    ES = write_solution(ES, m, "nonconvex")

    return ES
end

function my_callback_function(cb_data, cb_where::Cint) #heuristic callback function Gurobi
    
    # callbacks are used for a variety of things (we are only interested in the MIPNODE )
    if cb_where != GRB_CB_MIPSOL && cb_where != GRB_CB_MIPNODE
        return
    end

    # we want to check if the node of the B&B tree has been solved to optimality
    if cb_where == GRB_CB_MIPNODE
        resultP = Ref{Cint}()
        GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
        if resultP[] != GRB_OPTIMAL
            return  # Solution is something other than optimal.
        end
    end

    Gurobi.load_callback_variable_primal(cb_data, cb_where)
    # println("callback called")
    @info "callback called"
    
    y_val = callback_value.(cb_data, m[:y])
    # y_hat_val = callback_value(cb_data, y_hat)
    xp_val = callback_value.(cb_data, m[:xp])
    xn_val =  callback_value.(cb_data, m[:xn])
    b_val = callback_value.(cb_data, m[:b])
    # println("y_val: ", y_val)

    #filter positive ReLUs and negative ReLUs that are not violated
    ReLU_mask_pos = y_val - xp_val - xn_val .== 0
    # println(sum(ReLU_mask_pos))
    @info "ReLU_mask_pos: $(sum(ReLU_mask_pos))"
    ReLU_mask_neg = xp_val .== 0
    @info "ReLU_mask_neg: $(sum(ReLU_mask_neg))"
    # println(sum(ReLU_mask_neg))

    #filter fractional binaries corresponding to these ReLUs
    b_mask_pos = 0 .< b_val[ReLU_mask_pos] .< 1
    # println(sum(b_mask_pos))
    @info "b_mask_pos: $(sum(b_mask_pos))"
    b_mask_neg = 0 .< b_val[ReLU_mask_neg] .< 1
    # println(sum(b_mask_neg))
    @info "b_mask_neg: $(sum(b_mask_neg))"

    #get the indices of the fractional binaries
    #get the variables corresponding to these indices
    b_frac_pos = m[:b][ReLU_mask_pos][b_mask_pos]
    b_frac_neg = m[:b][ReLU_mask_neg][b_mask_neg]

    @info "Fractional ReLUs: positive $(sum(b_mask_pos)), negative $(sum(b_mask_neg)), both $(sum(b_mask_both)) "
    b_var_sel = [b_frac_pos; b_frac_neg]
    # b_var_val = [floor.(Int,b_val[ReLU_mask_pos][b_mask_pos]); ceil.(Int,b_val[ReLU_mask_neg][b_mask_neg])] #might need to be float instead of int
    b_var_val = [floor.(b_val[ReLU_mask_pos][b_mask_pos]); ceil.(Int,b_val[ReLU_mask_neg][b_mask_neg])] #might need to be float instead of int

    # println("b_var_sel: ", b_var_sel)
    # println("b_var_val: ", b_var_val)
    @info "b_var_sel: $(b_var_sel)"
    @info "b_var_val: $(b_var_val)"
    if isempty(b_var_sel)
        # println("No fractional ReLUs")
        # println("Skipping")
        @info "No fractional ReLUs"
        @info "Skipping"
        return
    else
        status = MOI.submit(
            m, MOI.HeuristicSolution(cb_data), b_var_sel, b_var_val
        ) #test if we need a vector or vector of vectors...
        # println("I submitted a heuristic solution, and the status was: ", status)
        @info "I submitted a heuristic solution, and the status was: $(status)"
    end

end




# function my_callback_function_heur(cb_data::CPLEX.CallbackContext, context_id::Clong)
#     if context_id == CPX_CALLBACKCONTEXT_RELAXATION
#         CPLEX.load_callback_variable_primal(cb_data, context_id)

#         y_val = callback_value.(cb_data, m[:y])
#         # y_hat_val = callback_value.(cb_data, m[:y_hat])
#         xp_val = callback_value.(cb_data, m[:xp])
#         xn_val =  callback_value.(cb_data, m[:xn])
#         b_val = callback_value.(cb_data, m[:b])
#         # ... heuristic stuff ...

#         ##########################################################
#         #--------- 1. Filter fractional binaries ----------------#
#         ##########################################################

#         #filter positive ReLUs and negative ReLUs that are not violated
#         ReLU_mask_pos = y_val - xp_val - xn_val .== 0
#         # println(sum(ReLU_mask_pos))
#         # @info "ReLU_mask_pos: $(sum(ReLU_mask_pos))"
#         ReLU_mask_neg = xp_val .== 0
#         # @info "ReLU_mask_neg: $(sum(ReLU_mask_neg))"
#         # println(sum(ReLU_mask_neg))

#         #filter fractional binaries corresponding to these ReLUs
#         b_mask_pos = 0 .< b_val[ReLU_mask_pos] .< 1
#         # println(sum(b_mask_pos))
#         # @info "b_mask_pos: $(sum(b_mask_pos))"
#         b_mask_neg = 0 .< b_val[ReLU_mask_neg] .< 1
#         # println(sum(b_mask_neg))
#         # @info "b_mask_neg: $(sum(b_mask_neg))"

#         #get the indices of the fractional binaries
#         #get the variables corresponding to these indices
#         b_frac_pos = m[:b][ReLU_mask_pos][b_mask_pos]
#         b_frac_neg = m[:b][ReLU_mask_neg][b_mask_neg]

#         @info "Fractional ReLUs: positive $(sum(b_mask_pos)) and negative $(sum(b_mask_neg)) "
#         b_var_sel = [b_frac_pos; b_frac_neg]
#         # b_var_val = [floor.(Int,b_val[ReLU_mask_pos][b_mask_pos]); ceil.(Int,b_val[ReLU_mask_neg][b_mask_neg])] #might need to be float instead of int
#         b_var_val = [floor.(b_val[ReLU_mask_pos][b_mask_pos]); ceil.(Int,b_val[ReLU_mask_neg][b_mask_neg])] #might need to be float instead of int

#         @info "b_var_sel: $(b_var_sel)"
#         @info "b_var_val: $(b_var_val)"
#         if isempty(b_var_sel)
#             # println("No fractional ReLUs")
#             # println("Skipping")
#             @info "No fractional ReLUs"
#             @info "Skipping"
#             return
#         else
#             # status = MOI.submit(
#             #     m, MOI.HeuristicSolution(cb_data), b_var_sel, b_var_val
#             # ) #test if we need a vector or vector of vectors...
#             @info "trying to set heuristic solution..."
#             status_heur = CPXcallbackpostheursoln(
#                 cb_data,
#                 Cint(length(b_var_sel)),
#                 Cint[CPLEX.column(cb_data, index(var)) - 1 for var in b_var_sel],
#                 b_var_val,
#                 503,
#                 CPXCALLBACKSOLUTION_PROPAGATE,
#             )
#             #CPXCALLBACKSOLUTION_PROPAGATE CPXCALLBACKSOLUTION_SOLVE
#             #PROPOGATE works better than SOLVE for some reason...   
#             #Cint[CPLEX.column(cb_data, index(var)) - 1 for var in b_var_sel]
#             #Cint[_info(model, var).column - 1 for var in b_var_sel]
#             #Cint[CPLEX._info(model, var).column - 1 for var in b_var_sel]
#             #CPLEX.column(cb_data, index(m[:b][b_chs]))
#             # println("I submitted a heuristic solution, and the status was: ", status)
#             @info "I submitted a heuristic solution, and the status was: $(status_heur)"
#         end
#         #CPXcallbackpostheursoln(context, cnt, ind, val, obj, strat)

#     elseif context_id == CPX_CALLBACKCONTEXT_BRANCHING
#         # ... branching stuff ...
#         CPLEX.load_callback_variable_primal(cb_data, context_id)

#         y_val = callback_value.(cb_data, m[:y])
#         y_hat_val = callback_value.(cb_data, m[:y_hat])
#         xp_val = callback_value.(cb_data, m[:xp])
#         xn_val =  callback_value.(cb_data, m[:xn])
#         b_val = callback_value.(cb_data, m[:b])

#         ##########################################################
#         #--------- 1. Filter fractional binaries ----------------#
#         ##########################################################

#         #filter positive ReLUs and negative ReLUs that are not violated
#         ReLU_mask_pos = y_val - xp_val - xn_val .== 0
#         # println(sum(ReLU_mask_pos))
#         # @info "ReLU_mask_pos: $(sum(ReLU_mask_pos))"
#         ReLU_mask_neg = xp_val .== 0
#         # @info "ReLU_mask_neg: $(sum(ReLU_mask_neg))"
#         # println(sum(ReLU_mask_neg))

#         #filter fractional binaries corresponding to these ReLUs
#         b_mask_pos = 0 .< b_val[ReLU_mask_pos] .< 1
#         # println(sum(b_mask_pos))
#         # @info "b_mask_pos: $(sum(b_mask_pos))"
#         b_mask_neg = 0 .< b_val[ReLU_mask_neg] .< 1
#         # println(sum(b_mask_neg))
#         # @info "b_mask_neg: $(sum(b_mask_neg))"

#         #get the indices of the fractional binaries
#         #get the variables corresponding to these indices
#         b_frac_pos = m[:b][ReLU_mask_pos][b_mask_pos]
#         b_frac_neg = m[:b][ReLU_mask_neg][b_mask_neg]
#         @info "\n\n testing to see if heuristic work \n\n"
#         @info "Fractional ReLUs: positive $(sum(b_mask_pos)) and negative $(sum(b_mask_neg)) "

        

#         ##########################################################
#         #--------- 2. Compute branching strategy ----------------#
#         ##########################################################

#         # d = (.-xn_val.+xp_val.+1)./(y_val.+1) #add the '+1' to avoid division by zero #old version
#         y_hat_val_p = min.(y_hat_val, 0)
#         d = (.-xn_val.+xp_val.+1)./(y_val-y_hat_val_p.+1)
#         # @info "d: $(d)"
#         b_chs = argmax(d)
#         # @info "b_chs: $(b_chs)"

#         column = CPLEX.column(cb_data, index(m[:b][b_chs]))-1
#         # @info "column: $(column)"
#         # column = CPLEX.column(cb_data, index(m[:x]))
#         seqnum_p_L = Ref{Cint}(0)
#         seqnum_p_U = Ref{Cint}(1)  

#         #create the child node where the variable is fixed to 1 (set lower bound to 1)
#         status_L = CPXcallbackmakebranch(
#             cb_data, 
#             1,              # varcnt
#             Cint[column],   # varind
#             Cchar['L'],     # varlu
#             Cdouble[1.0],   # varbd #need to make two child nodes... 
#             0,              # rcnt
#             0,              # nzcnt
#             Cdouble[],      # rhs
#             Cchar[],        # sense
#             Cint[],         # rmatbeg
#             Cint[],         # rmatind
#             Cdouble[],      # rmatval
#             503,            # nodeest  TODO: what is an objective estimate for the branch?
#             seqnum_p_L,       # seqnum_p
#         )
#         #create the child node where the variable is fixed to 0 (set upper bound to 0)
#         status_U = CPXcallbackmakebranch(
#             cb_data, 
#             1,              # varcnt
#             Cint[column],   # varind
#             Cchar['U'],     # varlu
#             Cdouble[0.0],   # varbd #need to make two child nodes... 
#             0,              # rcnt
#             0,              # nzcnt
#             Cdouble[],      # rhs
#             Cchar[],        # sense
#             Cint[],         # rmatbeg
#             Cint[],         # rmatind
#             Cdouble[],      # rmatval
#             503,            # nodeest  TODO: what is an objective estimate for the branch?
#             seqnum_p_U,       # seqnum_p
#         )

#         if status_L != 0
#             @warn "CPXcallbackmakebranch failed with status $(status_L)"
#         elseif status_U != 0
#             @warn "CPXcallbackmakebranch failed with status $(status_U)"
#         else
#             @info "$(column)"
#             @info "Created two new branch that branches on $(b_chs).\n The status_U was: $(status_U) (successful). \n The status_L was: $(status_L) (successful)"
        
#         end
#     else
#         # ... other stuff ...
#         return
#     end
# end

function my_callback_function_vf(cb_data::CPLEX.CallbackContext, context_id::Clong)# CPLEX callback + variable fixing
    if context_id != CPX_CALLBACKCONTEXT_BRANCHING
      return
    end
    # Let's assume I cant to set the upper bound of m[:x] to 1.0
    # @info "callback called"
    CPLEX.load_callback_variable_primal(cb_data, context_id)

    y_val = callback_value.(cb_data, m[:y])
    y_hat_val = callback_value.(cb_data, m[:y_hat])
    xp_val = callback_value.(cb_data, m[:xp])
    xn_val =  callback_value.(cb_data, m[:xn])
    b_val = callback_value.(cb_data, m[:b])

    ##########################################################
    #--------- 1. Filter fractional binaries ----------------#
    ##########################################################

    #filter positive ReLUs and negative ReLUs that are not violated
    # find ReLUs that are not violated where negative part is zero and positive part is not zero #Maybe do between +/- 1e-6 instead of zero?
    ReLU_mask_pos = xn_val .== 0 .&& xp_val .> 0 # if the negative part is zero, then the ReLU is not violated, as the positive part must be correct. Maybe do between +/- 1e-6 instead of zero?
    
    #Find ReLUs that are not violated where positive part is zero and negative part is not zero
    ReLU_mask_neg = xp_val .== 0 .&& xn_val .< 0
    # find ReLUs that are not violated and both positive and negative parts are zero
    ReLU_mask_both = xp_val .== xn_val .== 0
    

    #filter fractional binaries corresponding to these ReLUs
    b_mask_pos = 0 .< b_val[ReLU_mask_pos] .< 1
    # println(sum(b_mask_pos))
    # @info "b_mask_pos: $(sum(b_mask_pos))"
    b_mask_neg = 0 .< b_val[ReLU_mask_neg] .< 1
    # println(sum(b_mask_neg))
    # @info "b_mask_neg: $(sum(b_mask_neg))"
    b_mask_both = 0 .< b_val[ReLU_mask_both] .< 1


    #get the indices of the fractional binaries
    #get the variables corresponding to these indices
    b_frac_pos = m[:b][ReLU_mask_pos][b_mask_pos]
    b_frac_neg = m[:b][ReLU_mask_neg][b_mask_neg]
    b_frac_both = m[:b][ReLU_mask_both][b_mask_both]
    #MAKE SURE THAT INDEX DOES NOT APPEAR TWICE, WHEN WE HAVE BOTH POSITIVE AND NEGATIVE RELU PARTS = 0 
    @info "Fractional ReLUs: positive $(sum(b_mask_pos)), negative $(sum(b_mask_neg)) and both $(sum(b_mask_both)) "
    # b_var_sel = [b_frac_pos; b_frac_neg] #
    b_var_sel = [b_frac_pos; b_frac_neg; b_frac_both]
    # b_var_val = [floor.(Int,b_val[ReLU_mask_pos][b_mask_pos]); ceil.(Int,b_val[ReLU_mask_neg][b_mask_neg])] #might need to be float instead of int
    b_var_val = [b_val[ReLU_mask_pos][b_mask_pos]; b_val[ReLU_mask_neg][b_mask_neg];b_val[ReLU_mask_both][b_mask_both]]
    
    # @info "B_both: $(b_frac_both)"
    # @info "\nPOSITIVE"
    # @info "Calculation for violated ReLU y_hat_val $(y_hat_val[ReLU_mask_pos][b_mask_pos])"
    # @info "Calculation for violated ReLu y_val: $(y_val[ReLU_mask_pos][b_mask_pos])"
    # @info "Calculation for violated ReLu xp_val: $(xp_val[ReLU_mask_pos][b_mask_pos])"
    # @info "Calculation for violated ReLu xn_val: $(xn_val[ReLU_mask_pos][b_mask_pos])"

    # @info "\nNEGATIVE"
    # @info "Calculation for violated ReLU y_hat_val $(y_hat_val[ReLU_mask_neg][b_mask_neg])"
    # @info "Calculation for violated ReLu y_val: $(y_val[ReLU_mask_neg][b_mask_neg])"
    # @info "Calculation for violated ReLu xp_val: $(xp_val[ReLU_mask_neg][b_mask_neg])"
    # @info "Calculation for violated ReLu xn_val: $(xn_val[ReLU_mask_neg][b_mask_neg])"

    # @info "\nBOTH"
    # @info "Calculation for violated ReLU y_hat_val $(y_hat_val[ReLU_mask_both][b_mask_both])"
    # @info "Calculation for violated ReLu y_val: $(y_val[ReLU_mask_both][b_mask_both])"
    # @info "Calculation for violated ReLu xp_val: $(xp_val[ReLU_mask_both][b_mask_both])"
    # @info "Calculation for violated ReLu xn_val: $(xn_val[ReLU_mask_both][b_mask_both])"

    @info "b_var_val before rounding: $(b_var_val)"
    # b_var_val = [floor.(Int,b_val[ReLU_mask_pos][b_mask_pos]); ceil.(Int,b_val[ReLU_mask_neg][b_mask_neg])]
    b_var_val = [floor.(Int,b_val[ReLU_mask_pos][b_mask_pos]); ceil.(Int,b_val[ReLU_mask_neg][b_mask_neg]); ceil.(Int, b_val[ReLU_mask_both][b_mask_both])] #might need to be float instead of int
    seqnum_p_F = Ref{Cint}(2)
    # @info "b_var_sel: $(b_var_sel)"
    # @info "b_var_val: $(b_var_val)"

    ##########################################################
    #--------- 2. Compute branching strategy ----------------#
    ##########################################################
    # d = (.-xn_val.+xp_val.+1)./(y_val.+1) #add the '+1' to avoid division by zero #old version
    y_hat_val_p = min.(y_hat_val, 0)
    # d = (.-xn_val.+xp_val.+1)./(y_val-y_hat_val_p.+1) #????
    b_frac = 0.0 .< b_val .< 1.0
    # b_frac = findall(0.0 .< b_val .< 1.0)

    # d = (.-xn_val[b_frac].+xp_val[b_frac].+1)./(y_val[b_frac]-y_hat_val_p[b_frac].+1)

    d_all = (.-xn_val.+xp_val.+1)./(y_val-y_hat_val_p.+1)

    d = d_all.*b_frac

    if !all(y->y<= 1.0 , d) #fix the fractional binaries in the current node in a seperate branch.
        @info "No fractional ReLUs to fix, using custom branching rule"     

        # @info "d: $(d[d .> 1.0])" Fractional ReLUs: positive $(sum(b_mask_pos)), negative $(sum(b_mask_neg)) and both $(sum(b_mask_both)) 
        @info "\n\n--- total fractional binaries: $(sum(b_frac))---\n---fractional binaries with wrong ReLU (from d): $(sum(d .> 1.0))--\n---fixed variables: $(sum(b_mask_pos)+sum(b_mask_neg)+sum(b_mask_both)) -> $(sum(b_frac)) - $(sum(b_mask_pos)+sum(b_mask_neg)+sum(b_mask_both)) = $(sum(b_frac) - (sum(b_mask_pos)+sum(b_mask_neg)+sum(b_mask_both))) --- \n less than 1.0: $(sum(0.5 .< d .< 1.0))"
        
        @info "d[0.5 .< d .< 1.0]: $(d[0.5 .< d .< 1.0])"
        @info "y_val[0.5 .< d .< 1.0]: $(y_val[0.5 .< d .< 1.0])"
        @info "y_hat_val[d .< 1.0]: $(y_hat_val[0.5 .< d .< 1.0])"
        @info "b_val[d .< 1.0]: $(b_val[0.5 .< d .< 1.0])\n\n"

        
        # if all(y->y<= 1.0 , d)
        #     @info "\n-----------------------------------------------\n all d are 1.0 \n-----------------------------------------------\n"
        #     return
        # end
        # @info "d: $(d)"
        b_chs = argmax(d)
        @info "\n\n---VARIABLE TO BRANCH ON:\nb_chs: $(b_chs), $(maximum(d))"
        @info "Values of y_hat_val$(b_chs): $(y_hat_val[b_chs])"
        @info "Values of y_val$(b_chs): $(y_val[b_chs])"
        @info "Values of xp_val$(b_chs): $(xp_val[b_chs])"
        @info "Values of xn_val$(b_chs): $(xn_val[b_chs])"
        @info "Values of b_val$(b_chs): $(b_val[b_chs])\n\n"

        # @info "\nREFERENCE INDEX: CartesianIndex(23, 1, 1, 4): \n b_val: $(b_val[CartesianIndex(23, 1, 1, 4)]), \n yval $(y_val[CartesianIndex(23, 1, 1, 4)]),\n y_hat_val $(y_hat_val[CartesianIndex(23, 1, 1, 4)]),\n xp_val $(xp_val[CartesianIndex(23, 1, 1, 4)]),\n xn_val $(xn_val[CartesianIndex(23, 1, 1, 4)])\n"
        
        column = CPLEX.column(cb_data, index(m[:b][b_chs]))-1 #this is the variable to branch on... 

        @info "branching on $(column)"

        b_chs = argmax(d)
        # @info "b_chs: $(b_chs)"

        column = CPLEX.column(cb_data, index(m[:b][b_chs]))-1 #this is the variable to branch on... 
        # @info "column: $(column)"
        # column = CPLEX.column(cb_data, index(m[:x]))
        seqnum_p_L = Ref{Cint}(0)
        seqnum_p_U = Ref{Cint}(1)  

        
        # fix the fractional binaries with correct ReLU as part of the branching
        #maybe use this syntax for something...
        #Cint[CPLEX.column(cb_data, index(var)) - 1 for var in b_var_sel]
        #find a way to put 'L' and 'U' in the correct spots... needs to match the b_var_val (0 or 1)

        # define variables before branch call...

        #create the child node where the variable is fixed to 1 (set lower bound to 1)
        status_L = CPXcallbackmakebranch(
            cb_data, 
            1,              # varcnt
            Cint[column],   # varind
            Cchar['L'],     # varlu
            Cdouble[1.0],   # varbd #need to make two child nodes... 
            0,              # rcnt
            0,              # nzcnt
            Cdouble[],      # rhs
            Cchar[],        # sense
            Cint[],         # rmatbeg
            Cint[],         # rmatind
            Cdouble[],      # rmatval
            503,            # nodeest  TODO: what is an objective estimate for the branch?
            seqnum_p_L,       # seqnum_p
        )
        #create the child node where the variable is fixed to 0 (set upper bound to 0)
        status_U = CPXcallbackmakebranch(
            cb_data, 
            1,              # varcnt
            Cint[column],   # varind
            Cchar['U'],     # varlu
            Cdouble[0.0],   # varbd #need to make two child nodes... 
            0,              # rcnt
            0,              # nzcnt
            Cdouble[],      # rhs
            Cchar[],        # sense
            Cint[],         # rmatbeg
            Cint[],         # rmatind
            Cdouble[],      # rmatval
            503,            # nodeest  TODO: what is an objective estimate for the branch?
            seqnum_p_U,       # seqnum_p
        )

        if status_L != 0
            @warn "CPXcallbackmakebranch failed with status $(status_L)"
        elseif status_U != 0
            @warn "CPXcallbackmakebranch failed with status $(status_U)"
        else
            @info "$(column)"
            @info "Created two new branch that branches on $(b_chs).\n The status_U was: $(status_U) (successful). \n The status_L was: $(status_L) (successful)"
            @info "-------------------------------------------------"
        end
    else
        @info "Fractional ReLUs: positive $(sum(b_mask_pos)), negative $(sum(b_mask_neg)) and both $(sum(b_mask_both)) "
        status_F = CPXcallbackmakebranch(
            cb_data, 
            Cint(length(b_var_sel)),              # varcnt
            Cint[CPLEX.column(cb_data, index(var)) - 1 for var in b_var_sel],   # varind
            Cchar[val == 1.0 ? 'L' : 'U' for val in b_var_val],     # varlu #write 'L' or 'U' depending on the value of the variable (1.0 or 0.0)
            Cdouble[val for val in b_var_val],   # varbd #need to make two child nodes... 
            0,              # rcnt
            0,              # nzcnt
            Cdouble[],      # rhs
            Cchar[],        # sense
            Cint[],         # rmatbeg
            Cint[],         # rmatind
            Cdouble[],      # rmatval
            503,            # nodeest  TODO: what is an objective estimate for the branch?
            seqnum_p_F,       # seqnum_p
        )
        if status_F != 0
            @warn "CPXcallbackmakebranch failed with status $(status_F)"
        else
            # @info "$(column)"
            @info "Created a single branch that fixes the variables $(b_var_sel).\n The status_F was: $(status_F) (successful).\n The seqnum_p_F was: $(seqnum_p_F[])"
            @info "---------------------------------------------------------\n\n"
        end
    
    end
    
    
    
    # d = ones(3,3)

    # b_sel = [m[:b][b_chs]; b_var_sel]
    # b_val_U = [0.0; b_var_val]
    # b_val_L = [1.0; b_var_val]

    # @info "\n b_var_sel: $(b_sel)"
    # @info "\n b_var_val_U: $(b_val_U)"
    # @info "\n b_var_val_L: $(b_val_L)"

    # @info "column: $(column)"
    # column = CPLEX.column(cb_data, index(m[:x]))
    # seqnum_p_L = Ref{Cint}(0)
    # seqnum_p_U = Ref{Cint}(1)  
    
    # @info "Cchar val: $(Cchar[val == 1.0 ? 'L' : 'U' for val in b_val_L])"

    # 1.0 -> 'L' and 0.0 -> 'U'
    #create the child node where the variable is fixed to 1 (set lower bound to 1)
    # status_L = CPXcallbackmakebranch(
    #     cb_data, 
    #     Cint(length(b_sel)),              # varcnt
    #     Cint[CPLEX.column(cb_data, index(var)) - 1 for var in b_sel],   # varind
    #     Cchar[val == 1.0 ? 'L' : 'U' for val in b_val_L],     # varlu #write 'L' or 'U' depending on the value of the variable (1.0 or 0.0)
    #     Cdouble[val for val in b_val_L],   # varbd #need to make two child nodes... 
    #     0,              # rcnt
    #     0,              # nzcnt
    #     Cdouble[],      # rhs
    #     Cchar[],        # sense
    #     Cint[],         # rmatbeg
    #     Cint[],         # rmatind
    #     Cdouble[],      # rmatval
    #     503,            # nodeest  TODO: what is an objective estimate for the branch?
    #     seqnum_p_L,       # seqnum_p
    # )
    # #create the child node where the variable is fixed to 0 (set upper bound to 0)
    # status_U = CPXcallbackmakebranch(
    #     cb_data, 
    #     Cint(length(b_sel)),              # varcnt
    #     Cint[CPLEX.column(cb_data, index(var)) - 1 for var in b_sel],   # varind
    #     Cchar[val == 1.0 ? 'L' : 'U' for val in b_val_U],     # varlu #write 'L' or 'U' depending on the value of the variable (1.0 or 0.0)
    #     Cdouble[val for val in b_val_U],   # varbd #need to make two child nodes... 
    #     0,              # rcnt
    #     0,              # nzcnt
    #     Cdouble[],      # rhs
    #     Cchar[],        # sense
    #     Cint[],         # rmatbeg
    #     Cint[],         # rmatind
    #     Cdouble[],      # rmatval
    #     503,            # nodeest  TODO: what is an objective estimate for the branch?
    #     seqnum_p_U,       # seqnum_p
    # )

    # if status_L != 0
    #     @warn "CPXcallbackmakebranch failed with status $(status_L)"
    # elseif status_U != 0
    #     @warn "CPXcallbackmakebranch failed with status $(status_U)"
    # else
    #     @info "$(column)"
    #     @info "Created two new branch that branches on $(b_frac[b_chs]).\n The status_U was: $(status_U) (successful). \n The status_L was: $(status_L) (successful)\n The seqnum_p_L was: $(seqnum_p_L[])\n The seqnum_p_U was: $(seqnum_p_U[])"
    #     @info "---------------------------------------------------------\n\n"
    # end
end


function my_callback_function(cb_data::CPLEX.CallbackContext, context_id::Clong) #CPLEX CALLBACK
    if context_id != CPX_CALLBACKCONTEXT_BRANCHING
      return
    end
    # Let's assume I cant to set the upper bound of m[:x] to 1.0
    # @info "callback called"
    CPLEX.load_callback_variable_primal(cb_data, context_id)

    y_val = callback_value.(cb_data, m[:y])
    y_hat_val = callback_value.(cb_data, m[:y_hat])
    xp_val = callback_value.(cb_data, m[:xp])
    xn_val =  callback_value.(cb_data, m[:xn])
    b_val = callback_value.(cb_data, m[:b])

    ##########################################################
    #--------- 1. Filter fractional binaries ----------------#
    ##########################################################

    #filter positive ReLUs and negative ReLUs that are not violated
    ReLU_mask_pos = y_val - xp_val - xn_val .== 0
    # println(sum(ReLU_mask_pos))
    # @info "ReLU_mask_pos: $(sum(ReLU_mask_pos))"
    ReLU_mask_neg = xp_val .== 0
    # @info "ReLU_mask_neg: $(sum(ReLU_mask_neg))"
    # println(sum(ReLU_mask_neg))

    #filter fractional binaries corresponding to these ReLUs
    b_mask_pos = 0 .< b_val[ReLU_mask_pos] .< 1
    # println(sum(b_mask_pos))
    # @info "b_mask_pos: $(sum(b_mask_pos))"
    b_mask_neg = 0 .< b_val[ReLU_mask_neg] .< 1
    # println(sum(b_mask_neg))
    # @info "b_mask_neg: $(sum(b_mask_neg))"

    #get the indices of the fractional binaries
    #get the variables corresponding to these indices
    b_frac_pos = m[:b][ReLU_mask_pos][b_mask_pos]
    b_frac_neg = m[:b][ReLU_mask_neg][b_mask_neg]

    # @info "Fractional ReLUs: positive $(sum(b_mask_pos)) and negative $(sum(b_mask_neg)) "
    b_var_sel = [b_frac_pos; b_frac_neg]
    # b_var_val = [floor.(Int,b_val[ReLU_mask_pos][b_mask_pos]); ceil.(Int,b_val[ReLU_mask_neg][b_mask_neg])] #might need to be float instead of int
    b_var_val = [floor.(b_val[ReLU_mask_pos][b_mask_pos]); ceil.(Int,b_val[ReLU_mask_neg][b_mask_neg])] #might need to be float instead of int

    # @info "b_var_sel: $(b_var_sel)"
    # @info "b_var_val: $(b_var_val)"
    b_frac = 0.0 .< b_val .< 1.0
    ##########################################################
    #--------- 2. Compute branching strategy ----------------#
    ##########################################################

    # y_val

    # d = (.-xn_val.+xp_val.+1)./(y_val.+1) #add the '+1' to avoid division by zero #old version
    y_hat_val_p = min.(y_hat_val, 0)
    d_all = (.-xn_val.+xp_val.+1)./(y_val-y_hat_val_p.+1)
    d = d_all.*b_frac #only consider the fractional binaries
    if all(y->y<= 1.0 , d)
        @info "\n-----------------------------------------------\n all d are 1.0 \n-----------------------------------------------\n"
        return
    end
    # @info "d: $(d)"
    b_chs = argmax(d)
    # @info "b_chs: $(b_chs)"

    column = CPLEX.column(cb_data, index(m[:b][b_chs]))-1 #this is the variable to branch on... 
    # @info "column: $(column)"
    # column = CPLEX.column(cb_data, index(m[:x]))
    seqnum_p_L = Ref{Cint}(0)
    seqnum_p_U = Ref{Cint}(1)  

    
    # fix the fractional binaries with correct ReLU as part of the branching
    #maybe use this syntax for something...
    #Cint[CPLEX.column(cb_data, index(var)) - 1 for var in b_var_sel]
    #find a way to put 'L' and 'U' in the correct spots... needs to match the b_var_val (0 or 1)

    # define variables before branch call...

    #create the child node where the variable is fixed to 1 (set lower bound to 1)
    status_L = CPXcallbackmakebranch(
        cb_data, 
        1,              # varcnt
        Cint[column],   # varind
        Cchar['L'],     # varlu
        Cdouble[1.0],   # varbd #need to make two child nodes... 
        0,              # rcnt
        0,              # nzcnt
        Cdouble[],      # rhs
        Cchar[],        # sense
        Cint[],         # rmatbeg
        Cint[],         # rmatind
        Cdouble[],      # rmatval
        503,            # nodeest  TODO: what is an objective estimate for the branch?
        seqnum_p_L,       # seqnum_p
    )
    #create the child node where the variable is fixed to 0 (set upper bound to 0)
    status_U = CPXcallbackmakebranch(
        cb_data, 
        1,              # varcnt
        Cint[column],   # varind
        Cchar['U'],     # varlu
        Cdouble[0.0],   # varbd #need to make two child nodes... 
        0,              # rcnt
        0,              # nzcnt
        Cdouble[],      # rhs
        Cchar[],        # sense
        Cint[],         # rmatbeg
        Cint[],         # rmatind
        Cdouble[],      # rmatval
        503,            # nodeest  TODO: what is an objective estimate for the branch?
        seqnum_p_U,       # seqnum_p
    )

    if status_L != 0
        @warn "CPXcallbackmakebranch failed with status $(status_L)"
    elseif status_U != 0
        @warn "CPXcallbackmakebranch failed with status $(status_U)"
    else
        @info "$(column)"
        @info "Created two new branch that branches on $(b_chs).\n The status_U was: $(status_U) (successful). \n The status_L was: $(status_L) (successful)"
        @info "-------------------------------------------------"
    end
end


function allinfo(a)
    type = typeof(a)
    display(fieldnames(type))
    methodswith(type)
 end

function set_logger(path)
    folder = joinpath(dirname(@__DIR__), "output", path)
    if isdir(folder) == false
        mkdir(folder)
    end
    file = joinpath(folder, "log.txt")
    logger = SimpleLogger(open(file, "w+"))
    global_logger(logger)
end 
#  [b[2,21,1,4], b[12,5,1,1], b[4,6,1,1], b[6,14,1,3], b[11,14,1,3], b[17,15,1,3], b[14,18,1,3], b[17,20,1,3], b[11,12,2,3], b[9,7,1,4], b[23,20,1,4], b[7,20,2,4], b[9,20,2,4], b[10,20,2,4], b[16,20,2,4], b[24,20,2,4], b[23,15,2,5], b[17,6,1,3], b[1,7,1,3], b[18,7,1,3], b[6,14,1,3], b[19,14,1,3], b[4,15,1,3], b[17,15,1,3], b[15,17,1,3], b[14,18,1,3], b[21,18,1,3], b[2,19,1,3], b[12,19,1,3], b[7,20,1,3], b[23,5,2,3], b[19,15,2,3], b[8,20,2,3], b[24,20,2,3], b[3,6,1,5], b[5,16,1,5]]

# (https://www.ibm.com/docs/en/cofz/12.10.0?topic=context-cpx-callbackcontext-branching)

# https://discourse.julialang.org/t/error-using-cpxcallbackgetrelaxationpoint-at-cplex-jl-solver-specifc-callback/86230


# function my_callback_function(cb_data::CPLEX.CallbackContext, context_id::Clong)
#     if context_id != CPX_CALLBACKCONTEXT_BRANCHING

#         @info "callback called"

#         # Before querying `callback_value`, you must call:
#         CPLEX.load_callback_variable_primal(cb_data, context_id)

#         # get the relevant values from the callback
#         y_val = callback_value.(cb_data, m[:y])
#         xp_val = callback_value.(cb_data, m[:xp])
#         xn_val =  callback_value.(cb_data, m[:xn])
#         b_val = callback_value.(cb_data, m[:b])

#         #Compute the maximum distance to the ReLU hyperplane
#         d = (.-xn_val.+xp_val.+1)./(y_val.+1) #add the '+1' to avoid division by zero
#         #get the index of the associated fractional binary (if the ReLU is not enforced it will be > 1 and otherwise = 1)
#         b_chs = argmax(d)

#         #load the variable i want to branch on with proper CPLEX syntax:
#         column = CPLEX.column(cb_data, index(m[:b][b_chs]))
  
#         status = CPXcallbackmakebranch(
#             cb_data, 
#             1,              # varcnt
#             Cint[column],   # varind
#             Cchar['U'],     # varlu
#             Cdouble[1.0],   # varbd #need to make two child nodes... 
#             0,              # rcnt
#             0,              # nzcnt
#             Cdouble[],      # rhs
#             Cchar[],        # sense
#             Cint[],         # rmatbeg
#             Cint[],         # rmatind
#             Cdouble[],      # rmatval
#             503,            # nodeest  TODO: what is an objective estimate for the branch? -> i set the objective to 503 as that is slightly above the objective of the root node (and i know the solution is around there)
#             seqnum_p,       # seqnum_p
#         )
#     if status != 0
#         @warn "CPXcallbackmakebranch failed with status $(status)"
#     else
#         @info "I submitted a new branching strategy that branches on $(b_chs), and the status was: $(status) (successful)"
#     end
# end