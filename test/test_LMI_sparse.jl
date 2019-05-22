using JuMP, LinearAlgebra, SparseArrays
using SCS
using MosekTools
const MOI = JuMP.MathOptInterface

# formulating LMI constraint with sparse matrices
#   https://github.com/JuliaOpt/JuMP.jl/issues/1950
# We are trying to solve:
# min c' x
# s.t.  A0 + A1 * x[1] + A2 * x[2] ⪰ 0

m = 2
A0 = sparse([1. 0; 0 3])
A1 = sparse([0. 1; 1 -2])
A2 = sparse([5. -1; -1 0])
B = [A1, A2]
c = [1.; 1]

# A0 = [1. 0; 0 3]
# A1 = [0. 1; 1 -2]
# A2 = [5. -1; -1 0]
# B = [A1, A2]
# c = [1.; 1]

# model = JuMP.Model(with_optimizer(SCS.Optimizer, verbose = 0));
model = JuMP.Model(with_optimizer(Mosek.Optimizer, QUIET=false))

@variable(model, x[1:m]);
@objective(model, Min, c' * x);

# formulate LMI constraint A0 + A1 x[1] + A2 x[2] in PSDCone
# works
K = A0 + sum(B[k] .* x[k] for k in 1:m)
@SDconstraint(model, con1,  Symmetric(K) >= 0);
# @constraint(model, con1,  Symmetric(K) in JuMP.PSDCone());

# throws error
# @constraint(model, con2,  A0 + sum(B[k] .* x[k] for k in 1:m) in JuMP.PSDCone());
# @SDconstraint(model, con1,  K >= 0);

JuMP.optimize!(model)
@show JuMP.objective_value(model)

@show JuMP.termination_status(model) == MOI.OPTIMAL
@show JuMP.primal_status(model) == MOI.FEASIBLE_POINT

L = JuMP.value.(K)
# L = A0 + sum(B[k] .* obj_x[k] for k in 1:m)
@show isposdef(L)
