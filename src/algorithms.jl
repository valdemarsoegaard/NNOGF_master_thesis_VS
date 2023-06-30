include("opt_models_gas.jl")
include("postprocessing.jl")


function SLP(ES::EnergySystem, k_max::Int)
    # Successive linear programming algorithms
    """
    
    Attribute
    k::Iteration number
    """

    Q_nm_k = 0 #zeros(length(ES.T), length(ES.P), k_max)
    pr_avg_k = 0# zeros(length(ES.T), length(ES.P), k_max) 

    objective_values = []

    m = build_PDE_model_SLP(ES)
    # Define objective function
    obj_func = @expression(m, ES.dt/3600 *(
        sum(m[:Q_w][t,s]*ES.sources[s].c1 + m[:Q_w][t,s]^2*ES.sources[s].c2 for s in ES.S, t in ES.T) +
        sum(m[:Qd_cur_gas][t,d]*ES.C_cur_load_gas for d in ES.D, t in ES.T) +
        # sum(0.001*m[:pr][t,n] for n in ES.N, t in ES.T) +
        (isempty(ES.compressors) ? 0.0 : sum(ES.compressors[c].compr_cost*(
            m[:pr][t,ES.compressors[c].stop]-m[:pr][t,ES.compressors[c].start]) for c in ES.C, t in ES.T))
        ))

    if typeof(ES) == IntegratedEnergySystem
        nonGFPPs = [g for g in ES.G if ES.generators[g].type == "non-NGFPP"]
        obj_func = @expression(m, obj_func + 
            ES.dt/3600 *(sum(ES.generators[g].c_el*m[:p][g,t] for g in nonGFPPs, t in ES.T)))
    end

    # Define SLP parameters
    δ = [0, 10^-4] # Penalty factor augmentation
    δ_update = 2 # Update factor for penalty
    δ_max = 10^5 # Maximum penalty parameter
    ϕ = 0 # Difference between Taylor-series approximation of last iteration and current model value
    start = Dates.now()
    for k in 1:k_max

        println("Start iteration $k.")
        if k == 1
            # Penalize error term for first iteration
            set_objective_function(m, obj_func+sum(m[:γ][t,p]^2 for t in ES.T for p in ES.P))
        else
            # Add Taylor-series-based constraint for momentum conservation constraint
            m = add_taylor_series_expansion(ES, m, k, Q_nm_k, pr_avg_k)
            # Add augmentation terms for iteration k
            set_objective_function(m, obj_func + sum(
                δ[k]*(m[:pr_avg][t,p] - pr_avg_k[t,p])^2 + δ[k]*(m[:Q_nm][t,p] - Q_nm_k[t,p])^2
                for t in ES.T, p in ES.P))
        end
        # Solve model and extract required optimial values
        solve_model!(m)
        Q_nm = value.(m[:Q_nm])
        pr_avg = value.(m[:pr_avg])
        γ_k  = value.(m[:γ])#[value(m[:γ][t,p]) for t in ES.T, p in ES.P]
        append!(objective_values, objective_value(m))

        # Check solution convergence
        if k >= 2
            ϕ = abs.(γ_k .- Q_nm.*abs.(Q_nm)./pr_avg)
            if all(ϕ .<= 10^-8)

                ES = write_solution(ES, m, ES.model_type)
                # overwrite objective value with list including all iterations
                ES.obj_val = objective_values
                # ES.solve_time

                println("Process converged after $k out of $k_max iterations.")
                break
            end
        end

        # Save optimal values for next iteration
        Q_nm_k = Q_nm 
        pr_avg_k = pr_avg 

        println("End iteration $k.")
        if k == 1
            set_objective_function(m, obj_func)
        elseif k < k_max
            append!(δ, min(δ_update*δ[k], δ_max))
        else
            println("Process did NOT converge after $k_max iterations.
            You may want to try a higher value of k_max.")
            # include possible metric to capture gap here.
        end

    end
    stop = Dates.now()
    println("Solution time: $((stop-start)/Millisecond(1000)) seconds.")

    return ES

end


function add_taylor_series_expansion(ES, m, k, Q_nm_k, pr_avg_k)
    """
    Add description.
    """
    # If it's not the first iteration, the Taylor-series Expansion from
    # the last iteration has to be deleted first to avoid naming issues.
    if k > 2
        for t in ES.T for p in ES.P
                delete(m, m[:taylor_series_expansion][t,p])
        end end
        unregister(m, :taylor_series_expansion)
    end
    # Add Taylor-series expansion for current iteration
    @constraint(m, taylor_series_expansion[t in ES.T, p in ES.P],
        m[:γ][t,p] == (abs(Q_nm_k[t,p]) * Q_nm_k[t,p])/pr_avg_k[t,p] +
        (abs(Q_nm_k[t,p]) + sign(Q_nm_k[t,p])*Q_nm_k[t,p])/pr_avg_k[t,p] * (m[:Q_nm][t,p] - Q_nm_k[t,p]) -
        Q_nm_k[t,p]*abs(Q_nm_k[t,p])/pr_avg_k[t,p]^2 * (m[:pr_avg][t,p] - pr_avg_k[t,p])
    )
    return m
end