include("JumpModelToMatrix.jl")
include("DecompMatrix.jl")

module DW_ColGen

export GAPMIP, DWColGenEasy

using JuMP, GLPK
using LinearAlgebra, SparseArrays
using ..DecompMatrix, ..JumpModelToMatrix

struct ExtremePointInSol
    sub::Int64
    weight::Float64
    varVals::Vector{Float64}
end

EPSVAL = 0.00001

function setupSub(sub::JuMP.Model, Asub, bsub, sense, subvars, varLB, varUB, vecIsInt)
    (mAsub, n) = size(Asub)
    vecSubVarLB = zeros(Float64,n)
    vecSubVarUB = zeros(Float64,n)

    @variable(sub, xVars[1:n])
    # define upper and lower bounds as well as integrality constraints
    for j=1:n
         set_lower_bound(xVars[j], varLB[subvars[j]])
         set_upper_bound(xVars[j], varUB[subvars[j]])
         if vecIsInt[subvars[j]]
             set_integer(xVars[j])
         end
     end

    # objective is not here. We define once dual variables become known
    @objective(sub, Max, 0 )

    for i=1:mAsub
        if sense[i]==LEQ
            @constraint(sub, sum(Asub[i,j]*xVars[j] for j=1:n) <= bsub[i] )
        elseif sense[i]==GEQ
            @constraint(sub, sum(Asub[i,j]*xVars[j] for j=1:n) >= bsub[i] )
        else
            @constraint(sub, sum(Asub[i,j]*xVars[j] for j=1:n) == bsub[i] )
        end
    end
    return xVars
end

function setupMaster(master::JuMP.Model, A0, b0, nSub, sense, minimization=true)
    (mA0, nA0) = size(A0[1])
    K = 1

    # In this case we do not use a starting set of extreme points.
    # we just use a dummy starting column

    @variable(master, lambda[1:K] >= 0 )
    # we always solve the problem as a minimization problem, even if it
    # originally is a maximization problem. In case of maximization we
    # flip the sign and modify the output (such that it seems like we are
    # solving a maximization problem to the user)
    @objective(master, Min, sum(1000 * lambda[j] for j=1:K) )
    JuMPConstraintRef = JuMP.ConstraintRef
    consref = Vector{JuMPConstraintRef}(undef,mA0)
    for i=1:mA0
        if sense[i]==LEQ
            consref[i] = @constraint(master, sum(b0[i]*lambda[j] for j=1:K) <= b0[i] )
        elseif sense[i]==GEQ
            consref[i] = @constraint(master, sum(b0[i]*lambda[j] for j=1:K) >= b0[i] )
        else
            consref[i] = @constraint(master, sum(b0[i]*lambda[j] for j=1:K) == b0[i] )
        end
    end
    @constraint(master, convexityCons[k=1:nSub], sum(lambda[j] for j=1:K) == 1)

    return consref, convexityCons, lambda
end

# x is the new extreme point we wish to add to the master problem
function addColumnToMaster(master::JuMP.Model, pPerSub, A0, x, consRef, convexityCons, subNumber, minimize::Bool)
    (mA0, n) = size(A0[subNumber])
    A0x = A0[subNumber]*x
    if minimize
        mult = 1
    else
        mult = -1
    end
    oldvars = JuMP.all_variables(master)
    new_var = @variable(master, base_name="lambda_$(length(oldvars))", lower_bound=0)
    JuMP.set_objective_coefficient(master, new_var, mult*sum(pPerSub[subNumber][j]*x[j] for j=1:n))

    for i=1:mA0
        # only insert non-zero elements (this saves memory and may make the master problem easier to solve)
        if A0x[i] != 0
            set_normalized_coefficient(consRef[i], new_var, A0x[i])
        end
    end
    # add variable to convexity constraint.
    set_normalized_coefficient(convexityCons[subNumber], new_var, 1)

    return new_var

end


function solveSub(sub, myPi, myKappa, pPerSub, A0, xVars,subNumber, minimize::Bool)
    pSub = pPerSub[subNumber]
    A0Sub = A0[subNumber]
    xVarsSub = xVars[subNumber]
    #println("typeof(A0): ", typeof(A0))
    #println("typeof(A0Sub): ", typeof(A0Sub))
    (mA0, n) = size(A0Sub)

    piA0 = myPi*A0Sub
    if minimize
        mult = 1
    else
        mult = -1
    end

    # set objective. Remember to consider if maximization or minimization is needed
    @objective(sub[subNumber], Min, mult*sum(pSub[j]*xVarsSub[j] for j=1:n) - sum(piA0[1,j]*xVarsSub[j] for j=1:n) - myKappa[subNumber])

    optimize!(sub[subNumber])
    if termination_status(sub[subNumber]) != MOI.OPTIMAL
        throw("Error: Non-optimal sub-problem status")
    end

    #println(sub[subNumber])
    #println("sub $subNumber objective: $(JuMP.objective_value(sub[subNumber])), vars: $(value.(xVarsSub))")

    return JuMP.objective_value(sub[subNumber]), value.(xVarsSub)
end

function DWColGen(A0,ASub,b0,bSub, senseA0, senseSub, pPerSub,subvars, mip::MIP)
    if mip.minimization
        mult = 1
    else
        mult = -1
    end

    varNames = mip.varNames
    varLB = mip.varLB
    varUB = mip.varUB
    vecIsInt = mip.vecIsInt

    nSub = length(bSub)
    sub = Vector{JuMP.Model}(undef, nSub)

    for k=1:nSub
        sub[k] = Model(GLPK.Optimizer)
    end
    master = Model(GLPK.Optimizer)
    xVars = []
    for k=1:nSub
        push!(xVars, setupSub(sub[k], ASub[k], bSub[k], senseSub[k], subvars[k], varLB, varUB, vecIsInt))
    end
    (consRef, convexityCons, lambdas) = setupMaster(master, A0, b0, nSub, senseA0, mip.minimization)
    # extremePoints records the extreme points. extremePoints[p] corresponds to variable lambda[p]
    # extremePointForSub[p] records which sub-problem the p'th extreme point "belongs to"
    # first extreme point is a dummy one.
    extremePoints = [[]]
    extremePointForSub = [-1]
    done = false
    iter = 1
    while !done
        optimize!(master)
        if termination_status(master) != MOI.OPTIMAL
            throw("Error: Non-optimal master-problem status")
        end
        myPi = dual.(consRef)
        # ensure that myPi and myKappa are  row vectors
        myPi = reshape(myPi, 1, length(myPi))
        myKappa = dual.(convexityCons)
        myKappa = reshape(myKappa, 1, length(myKappa))
        #println("myPi = $myPi")
        #println("myKappa = $myKappa")
        done = true
        println("iteration: $iter, objective value = $(mult*JuMP.objective_value(master))")
        bestRedCost = 1
        for k=1:nSub
            redCost, xVal = solveSub(sub, myPi, myKappa, pPerSub, A0, xVars, k, mip.minimization)
            print("sub $k red cost = $(mult*redCost), ")
            if redCost < bestRedCost
                bestRedCost = redCost
            end
            if redCost < -EPSVAL
                newVar = addColumnToMaster(master, pPerSub, A0, xVal, consRef, convexityCons,k,mip.minimization)
                push!(lambdas,newVar)
                push!(extremePoints, xVal)
                push!(extremePointForSub, k)
                done = false
            end
        end
        iter += 1
        println("best reduced cost = $(mult*bestRedCost)")
    end
    println("Done after $(iter-1) iterations. Objective value = $(mult*JuMP.objective_value(master))")
    # compute values of original variables
    origVarValSub = []
    for s = 1:length(subvars)
        push!(origVarValSub,zeros(length(subvars[s])))
    end
    lambdaVal = value.(lambdas)
    usedExtremePoints = Vector{ExtremePointInSol}()
    for p=1:length(lambdaVal)
        if lambdaVal[p] > 0.0001
            println("lambda_$p=", lambdaVal[p], ", sub=$(extremePointForSub[p]), extr.point=$(extremePoints[p])")
            origVarValSub[extremePointForSub[p]] += lambdaVal[p]*extremePoints[p]
            push!(usedExtremePoints, ExtremePointInSol(extremePointForSub[p], lambdaVal[p], extremePoints[p]))
        end
    end
    vecOrigVars = zeros(Float64, length(varNames))
    for s = 1:length(subvars)
        #println("var val for sub problem $s: $(origVarValSub[s])")
        for t=1:length(origVarValSub[s])
            if abs(origVarValSub[s][t]) > 0.0001
                vecOrigVars[subvars[s][t]] = origVarValSub[s][t]
                println("$(varNames[subvars[s][t]])=$(origVarValSub[s][t])")
            end
        end
    end
    return usedExtremePoints, vecOrigVars
end

function DWColGenEasy(mip::MIP, blocks::Vector{Vector{Int64}})
    timeStart = time()
    A0, b0, A0Sense, Asub, bSub, subSense, subVars, cPerSub = constructSubMatrices(mip,blocks)
    usedExtremePoints, vecOrigVars = DWColGen(A0,Asub,b0,bSub, A0Sense, subSense, cPerSub,subVars,mip)
    #println("usedExtremePoints: ", usedExtremePoints)
    println("Elapsed time: $(time()-timeStart) seconds")
    return usedExtremePoints, vecOrigVars
end

function DWColGenEasy(jumpModel, blocks::Vector{Vector{ConstraintRef}})
    mip, constraintRefToRowIdDict = getConstraintMatrix(jumpModel)
    blockInt = Vector{Vector{Int64}}()
    for block in blocks
        push!(blockInt, Vector{Int64}())
        for constraintRef in block
            push!(blockInt[end], constraintRefToRowIdDict[constraintRef])
        end
        println("block $(length(blockInt)): ", blockInt[end])
    end

    return DWColGenEasy(mip, blockInt)
end

end
