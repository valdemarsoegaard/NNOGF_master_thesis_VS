using Bonobo
using JuMP
using HiGHS

const BB = Bonobo


m = Model(HiGHS.Optimizer)
set_optimizer_attribute(m, "log_to_console", false)
@variable(m, x[1:3] >= 0)
@constraint(m, 0.5x[1]+3.1x[2]+4.2x[3] >= 6.1)
@constraint(m, 1.9x[1]+0.7x[2]+0.2x[3] >= 8.1)
@constraint(m, 2.9x[1]-2.3x[2]+4.2x[3] >= 10.5)
@objective(m, Min, x[1]+1.2x[2]+3.2x[3])



mutable struct MIPNode <: AbstractNode
    std :: BnBNodeInfo
    lbs :: Vector{Float64}
    ubs :: Vector{Float64}
    status :: MOI.TerminationStatusCode
end


function BB.evaluate_node!(tree::BnBTree{MIPNode, JuMP.Model}, node::MIPNode)
    m = tree.root # this is the JuMP.Model
    vids = MOI.get(m ,MOI.ListOfVariableIndices())
    # we set the bounds for the current node based on `node.lbs` and `node.ubs`.
    vars = VariableRef.(m, vids)
    for vidx in eachindex(vars)
        if isfinite(node.lbs[vidx])
            JuMP.set_lower_bound(vars[vidx], node.lbs[vidx])
        elseif node.lbs[vidx] == -Inf && JuMP.has_lower_bound(vars[vidx])
            JuMP.delete_lower_bound(vars[vidx])
        elseif node.lbs[vidx] == Inf # making problem infeasible
            error("Invalid lower bound for variable $vidx: $(node.lbs[vidx])")
        end
        if isfinite(node.ubs[vidx])
            JuMP.set_upper_bound(vars[vidx], node.ubs[vidx])
        elseif node.ubs[vidx] == Inf && JuMP.has_upper_bound(vars[vidx])
            JuMP.delete_upper_bound(vars[vidx])
        elseif node.ubs[vidx] == -Inf # making problem infeasible
            error("Invalid upper bound for variable $vidx: $(node.lbs[vidx])")
        end
    end

    # get the relaxed solution of the current model using HiGHS
    optimize!(m)
    status = termination_status(m)
    node.status = status
    # if it is infeasible we return `NaN` for bother lower and upper bound
    if status != MOI.OPTIMAL
        return NaN,NaN
    end

    obj_val = objective_value(m)
    # we check whether the values are approximately feasible (are integer)
    # in that case we return the same value for lower and upper bound for this node
    if all(BB.is_approx_feasible.(tree, value.(vars)))
        node.ub = obj_val
        return obj_val, obj_val
    end
    # otherwise we only have a lower bound
    return obj_val, NaN
end

function BB.get_branching_indices(model::JuMP.Model)
    # every variable should be discrete
    vis = MOI.get(model, MOI.ListOfVariableIndices())
    return 1:length(vis)
end

function BB.get_relaxed_values(tree::BnBTree{MIPNode, JuMP.Model}, node)
    vids = MOI.get(tree.root, MOI.ListOfVariableIndices())
    vars = VariableRef.(tree.root, vids)
    return JuMP.value.(vars)
end

function BB.get_branching_nodes_info(tree::BnBTree{MIPNode, JuMP.Model}, node::MIPNode, vidx::Int)
    m = tree.root
    node_info = NamedTuple[]

    var = VariableRef(m, MOI.VariableIndex(vidx))

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

Bonobo.initialize
Bonobo.MOST_INFEASIBLE()

bnb_model = BB.initialize(;
    branch_strategy = BB.MOST_INFEASIBLE(),
    Node = MIPNode,
    root = m,
    sense = objective_sense(m) == MOI.MAX_SENSE ? :Max : :Min
)




BB.set_root!(bnb_model, (
    lbs = zeros(length(x)),
    ubs = fill(Inf, length(x)),
    status = MOI.OPTIMIZE_NOT_CALLED
))

BB.optimize!(bnb_model)

BB.get_solution(bnb_model)

BB.get_objective_value(bnb_model)


# b_var_sel: VariableRef[[22,9,1,4], [18,8,2,2], [23,9,1,3], b[19,7,1,4], b[18,8,1,4], b[7,19,1,4], b[18,20,2,4], b[6,1,1,1], b[10,1,1,1], b[18,1,1,1], b[6,2,1,1], b[10,2,1,1], b[18,2,1,1], b[4,3,1,1], b[4,4,1,1], b[6,4,1,1], b[15,5,1,1], b[8,6,1,1], b[2,7,1,1], b[20,7,1,1], b[22,9,1,1], b[3,10,1,1], b[9,10,1,1], b[20,10,1,1], b[22,10,1,1], b[3,11,1,1], b[9,11,1,1], b[12,11,1,1], b[20,11,1,1], b[22,11,1,1], b[3,14,1,1], b[22,16,1,1], b[22,17,1,1], b[18,18,1,1], b[3,19,1,1], b[17,19,1,1], b[2,21,1,1], b[11,20,1,2], b[5,6,1,3], b[9,6,1,3], b[7,7,1,3], b[11,12,1,3], b[13,14,1,3], b[22,14,1,3], b[8,15,1,3], b[19,15,1,3], b[23,15,1,3], b[19,17,1,3], b[1,19,1,3], b[6,19,1,3], b[7,19,1,3], b[16,19,1,3], b[19,19,1,3], b[22,19,1,3], b[2,20,1,3], b[24,20,1,3], b[24,4,2,3], b[1,5,2,3], b[13,5,2,3], b[4,9,2,3], b[22,9,2,3], b[12,12,2,3], b[17,14,2,3], b[14,15,2,3], b[22,15,2,3], b[2,16,2,3], b[19,17,2,3], b[8,18,2,3], b[18,20,2,3], b[23,20,2,3], b[6,12,1,5], b[6,13,1,5], b[3,16,1,5], b[3,18,1,5], b[1,13,2,5], b[5,13,2,5], b[6,18,2,5], b[18,20,2,5], b[22,20,2,5]]


A = [[22,9,1,4], [18,8,2,2], [23,9,1,3], [19,7,1,4], [18,8,1,4], [7,19,1,4], [18,20,2,4], [6,1,1,1], [10,1,1,1], [18,1,1,1], [6,2,1,1], [10,2,1,1], [18,2,1,1], [4,3,1,1], [4,4,1,1], [6,4,1,1], [15,5,1,1], [8,6,1,1], [2,7,1,1], [20,7,1,1], [22,9,1,1], [3,10,1,1], [9,10,1,1], [20,10,1,1], [22,10,1,1], [3,11,1,1], [9,11,1,1], [12,11,1,1], [20,11,1,1], [22,11,1,1], [3,14,1,1], [22,16,1,1], [22,17,1,1], [18,18,1,1], [3,19,1,1], [17,19,1,1], [2,21,1,1], [11,20,1,2], [5,6,1,3], [9,6,1,3], [7,7,1,3], [11,12,1,3], [13,14,1,3], [22,14,1,3], [8,15,1,3], [19,15,1,3], [23,15,1,3], [19,17,1,3], [1,19,1,3], [6,19,1,3], [7,19,1,3], [16,19,1,3], [19,19,1,3], [22,19,1,3], [2,20,1,3], [24,20,1,3], [24,4,2,3], [1,5,2,3], [13,5,2,3], [4,9,2,3], [22,9,2,3], [12,12,2,3], [17,14,2,3], [14,15,2,3], [22,15,2,3], [2,16,2,3], [19,17,2,3], [8,18,2,3], [18,20,2,3], [23,20,2,3], [6,12,1,5], [6,13,1,5], [3,16,1,5], [3,18,1,5], [1,13,2,5], [5,13,2,5], [6,18,2,5], [18,20,2,5], [22,20,2,5]]
B = hcat([[i, count(==(i), A)] for i in unique(A)]...)


#check unique elements in the list above... 


