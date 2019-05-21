using StabTrussTopOpt
sto = StabTrussTopOpt
using JuMP, SCS

function TO_solve_relaxed_stab(pb::sto.TOProblem)
    # here we solve the following relaxed layout opt problem w/ global stability:
    # (single load case)
    #       min_{a, q} l' * a
    # s.t.  ∑_{i ∈ [m]} q_i γ_i = f  (force equilibrium)
    #       - a σ^- ≦ q ≦ a σ^+      (stress)
    #       K(a) + τ G(q) ⪰ 0        (global stability)
    #       a ≧ 0
    #
    # JuMP modeling reference:
    #   https://github.com/JuliaOpt/JuMP.jl/blob/master/examples/robust_uncertainty.jl

    # scalar parameter
    m = pb.m # num of element
    E = pb.E
    σ_c = pb.σ_c
    σ_t = pb.σ_t
    τ = pb.τ # stability load factor
    nfree = size(pb.eK[1], 1)

    # vector & spmat data
    L = pb.l
    Γ = pb.Γ # ∈ m x nfree

    # sparse matrices array indexed by element id
    eK = pb.eK
    Δ = pb.Δ

    model = JuMP.Model(with_optimizer(SCS.Optimizer, verbose = 1))
    @variable(model, a[1:m] >= 0)
    @variable(model, q[1:m])

    @constraint(model, Γ'*q == f)

    K = spzeros(nfree, nfree)
    for i=1:m
        K += a[i] * eK[i]
    end
    G = spzeros(nfree, nfree)
    for i=1:m
        G += q[i]/L[i] * reshape(Δ[:,i], nfree, nfree)
    end
    @SDconstraint(model, K + τ * G >= 0)

    @constraint(model, - σ_c * a - q <= 0)
    @constraint(model, - σ_t * a + q <= 0)

    @objective(model, Min, dot(pb.c, a))

    JuMP.optimize!(model)
    @show JuMP.objective_value(model)
end
