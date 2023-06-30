# #################### Calculating pressures ####################
function calculate_pressures_from_ϕ(ES)
    pr = zeros(maximum(T), maximum(N))
    for t in T, p in P
        pr[t, pipes[p].start] = 1/2*(ES.ϕ⁺[t,p]+ES.ϕ⁻[t,p])
        pr[t, pipes[p].stop] = 1/2*(ES.ϕ⁺[t,p]-ES.ϕ⁻[t,p])
    end
    return pr
end

# #################### Calculating Feasibility Gap ####################
# function calculate_feasibility_gap(ES)
#     gap = [ES.Q_nm[t,p]^2 - ES.pipes[p].Knm^2 * (ES.pr[t,ES.pipes[p].start]^2 - ES.pr[t,ES.pipes[p].stop]^2)
#         for t in ES.T, p in ES.P]
#     return gap
# end

function calculate_feasibility_gap(ES)
    # if ES.model_type == "NN_constrained"
    gap = [value(ES.model[:Q_nm][t,p])*abs(value(ES.model[:Q_nm][t,p])) - ES.pipes[p].Knm^2 * (value(ES.model[:pr][t,ES.pipes[p].start])^2 - value(ES.model[:pr][t,ES.pipes[p].stop])^2)
            for t in ES.T, p in ES.P] #use Q*|Q| for NN
    # else
    #     gap = [value(ES.model[:Q_nm][t,p])^2 - ES.pipes[p].Knm^2 * (value(ES.model[:pr][t,ES.pipes[p].start])^2 - value(ES.model[:pr][t,ES.pipes[p].stop])^2)
    #                 for t in ES.T, p in ES.P] #use Q^2 for nonconvex and other models
    # end
    return round.(gap, digits = 5)
end



# ES.model[:Q_nm][t,p])^2 -
# [value(ES.pipes[p].Knm^2 * (value(ES.model[:pr][t,ES.pipes[p].start])^2 - value(ES.model[:pr][t,ES.pipes[p].stop])^2) for t in ES.T, p in ES.P]

# [ ES.pipes[p].Knm^2 * (value(ES.model[:pr][t,ES.pipes[p].start])^2 - value(ES.model[:pr][t,ES.pipes[p].stop])^2)
#         for t in ES.T, p in ES.P]

# #compute squared difference
# [(value(ES.model[:pr][t,ES.pipes[p].start])^2 - value(ES.model[:pr][t,ES.pipes[p].stop])^2) for t in ES.T, p in ES.P]

# [ES.pipes[p].Knm^2*(value(ES.model[:pr][t,ES.pipes[p].start])^2 - value(ES.model[:pr][t,ES.pipes[p].stop])^2) for t in ES.T, p in ES.P]
# # Comparare with ES.Q_nm
# [sqrt(ES.pipes[p].Knm^2*(value(ES.model[:pr][t,ES.pipes[p].start])^2 - value(ES.model[:pr][t,ES.pipes[p].stop])^2)) for t in ES.T, p in ES.P]