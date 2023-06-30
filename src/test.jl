using Pkg
using HiGHS
using JuMP
using Cbc
using Bonobo

# m = Model(HiGHS.Optimizer)
m = Model(Cbc.Optimizer)
set_optimizer_attribute(m, "logLevel", 0)
@variable(m, x[1:3] >= 0)
@constraint(m, 0.5x[1]+3.1x[2]+4.2x[3] <= 6.1)   
@constraint(m, 1.9x[1]+0.7x[2]+0.2x[3] <= 8.1)   
@constraint(m, 2.9x[1]-2.3x[2]+4.2x[3] <= 10.5)   
@objective(m, Max, x[1]+1.2x[2]+3.2x[3])

#functions needed

mutable struct MIPNode <: AbstractNode
    std :: BnBNodeInfo
    lbs :: Vector{Float64}
    ubs :: Vector{Float64}
    status :: MOI.TerminationStatusCode
end

function Bonobo.evaluate_node!(tree::BnBTree{MIPNode, JuMP.Model}, node::MIPNode)
    m = tree.root
    vids = MOI.get(m ,MOI.ListOfVariableIndices())
    vars = VariableRef.(m, vids)
    JuMP.set_lower_bound.(vars, node.lbs)
    JuMP.set_upper_bound.(vars, node.ubs)

    optimize!(m)
    status = termination_status(m)
    node.status = status
    if status != MOI.OPTIMAL
        return NaN,NaN
    end

    obj_val = objective_value(m)
    if all(Bonobo.is_approx_feasible.(tree, value.(vars)))
        node.ub = obj_val
        return obj_val, obj_val
    end
    return obj_val, NaN
end

function Bonobo.get_relaxed_values(tree::BnBTree{MIPNode, JuMP.Model}, node)
    vids = MOI.get(tree.root ,MOI.ListOfVariableIndices())
    vars = VariableRef.(tree.root, vids)
    return JuMP.value.(vars)
end

function Bonobo.get_branching_indices(model::JuMP.Model)
    # every variable should be discrete
    vis = MOI.get(model, MOI.ListOfVariableIndices())
    return 1:length(vis)
end

function Bonobo.get_branching_nodes_info(tree::BnBTree{MIPNode, JuMP.Model}, node::MIPNode, vidx::Int)
    m = tree.root
    node_info = NamedTuple[]

    var = VariableRef(m, MOI.VariableIndex(vidx))

    # first variable which is not discrete
    lbs = copy(node.lbs)
    ubs = copy(node.ubs)

    val = JuMP.value(var)

    # left child set upper bound
    ubs[vidx] = floor(Int, val)

    push!(node_info, (
        lbs = copy(node.lbs),
        ubs = ubs,
        status = MOI.OPTIMIZE_NOT_CALLED,
    ))

    # right child set lower bound
    lbs[vidx] = ceil(Int, val)

    push!(node_info, (
        lbs = lbs,
        ubs = copy(node.ubs),
        status = MOI.OPTIMIZE_NOT_CALLED,
    ))
    return node_info
end



bnb_model = Bonobo.initialize(; 
    Node = MIPNode,
    root = m,
    sense = objective_sense(m) == MOI.MAX_SENSE ? :Max : :Min
)

Bonobo.set_root!(bnb_model, (
    lbs = zeros(length(x)),
    ubs = fill(Inf, length(x)),
    status = MOI.OPTIMIZE_NOT_CALLED
))


Bonobo.optimize!(bnb_model)




Bonobo.BnBNode



2+2




#Psuedo code for branch and bound algorithm

#function branch_and_bound(m::Model)
#    # Solve initial LP relaxed model
#    JuMP.optimize!(m)
#    # Get the objective value of the relaxed model
#    obj_val = JuMP.objective_value(m)
#    # Get the value of the variables
#    x_val = JuMP.value.(x)
#    # Get the status of the relaxed model
#    status = JuMP.termination_status(m)
#    # If the relaxed model is infeasible, return infeasible
#    if status == MOI.INFEASIBLE
#        return "infeasible"
#    # If the relaxed model is feasible, return the objective value
#    elseif status == MOI.OPTIMAL
#        return obj_val
#    # If the relaxed model is not optimal, branch
#    else
        # Get the index of the variable with largest distance from the ReLU hyperplane
        # y_val = callback_value.(cb_data, m[:y])
        # # y_hat_val = callback_value(cb_data, y_hat)
        # xp_val = callback_value.(cb_data, m[:xp])
        # xn_val =  callback_value.(cb_data, m[:xn])
        # b_val = callback_value.(cb_data, m[:b])

        # ReLU_mask_pos = y_val - xp_val - xn_val .== 0

        # d = (-xn +xp+1)./(y+1)



#Use Bonobo.jl to solve the problem
