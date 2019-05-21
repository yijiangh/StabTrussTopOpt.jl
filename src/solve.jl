using StabTrussTopOpt
sto = StabTrussTopOpt

function TO_solve_relaxed_stab(par::sto.TOProblem)
    # here we solve the following relaxed layout opt problem w/ global stability:
    # (single load case)
    #       min_{a, q} l' * a
    # s.t.  ∑_{i ∈ [m]} q_i γ_i = f  (force equilibrium)
    #       - a σ^- ≦ q ≦ a σ^+      (stress)
    #       K(a) + τ G(q) ⪰ 0        (global stability)
    #       a ≧ 0

end
