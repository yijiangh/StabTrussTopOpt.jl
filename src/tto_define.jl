using StabTrussTopOpt
sto = StabTrussTopOpt

function TO_define(par::sto.TOProblemPar; E::Float64=1.0, τ::Float64=1.0, σ_c = 10.0, σ_t = 10.0)::sto.TOProblem

    ndim = par.ndim
    n = par.n
    m = par.m
    B = par.B
    edof_map = par.dof_map

    ndof = ndim * n
    # free dof
    mask_free = [v[1] for v in findall(i->i==0, reshape(par.Bound', ndof, 1))]
    nfree = length(mask_free)

    ff = reshape(par.F', ndof, 1)[mask_free]

    l = zeros(Float64, m)
    for i=1:m
        u_id = par.T[i,1]
        v_id = par.T[i,2]
        l[i] = norm(par.X[v_id,:] - par.X[u_id,:])
    end

    # linear objective coeff
    c = l

    g_I = Int[]
    g_J = Int[]
    g_V = Float64[]
    for i=1:nfree
        append!(g_I, i * ones(ndim))
        append!(g_J, edof_map[i, 1:ndim])
        append!(g_V, B[i, 1:ndim])

        append!(g_I, i * ones(ndim))
        append!(g_J, edof_map[i, ndim+1:2*ndim])
        append!(g_V, B[i, ndim+1:2*ndim])
    end
    Γ_full = sparse(g_I, g_J, g_V, m, ndof)
    Γ = Γ_full[:, mask_free] # ∈ m x nfree

    eK = Array{SparseMatrixCSC{Float64, Int}, 1}(undef, m)
    for i=1:m
        eK[i] = E / l[i] * (Γ[i, :] * Γ[i, :]')
    end

    # construct orthonormal basis for gamma_i, i ∈ [m]
    # aka. δ_i (and η_i in 3D)
    Bn = zeros(Float64, m, ndof)
    for i=1:m
        # TODO: general in 3D: Gram-Schmidt in the null space of gamma_i
        # here 2D specific construction to make orthonormal
        Bn[i,1] = B[i,2]
        Bn[i,2] = -B[i,1]
        Bn[i,3] = B[i,4]
        Bn[i,4] = -B[i,3]
    end

    # map Bn's local index 1..2xndim to global dof indices
    spBn_I = Int[]
    spBn_J = Int[]
    spBn_V = Float64[]
    for i=1:m
        append!(spBn_I, i * ones(ndim))
        append!(spBn_J, edof_map[i, 1:ndim])
        append!(spBn_V, Bn[i, 1:ndim])

        append!(spBn_I, i * ones(ndim))
        append!(spBn_J, edof_map[i, ndim+1:2*ndim])
        append!(spBn_V, Bn[i, ndim+1:2*ndim])
    end
    spBn = sparse(spBn_I, spBn_J, spBn_V, m, ndof)

    # construct Δ_i = (δ_i * δ_i^T + η_i * η_i^T), i ∈ [m]
    spDelta = spzeros(nfree^2, m)
    for i=1:m
        DDD = spBn[i, mask_free] * spBn[i, mask_free]'
        spDelta[:, i] = DDD[:]
    end

    a_init = 10 * ones(Float64, m)
    # TODO: solve the full comp-geometric stiffness matrix to get init_q
    q_init = 10 * ones(Float64, m)
    x_init  = vcat(a_init, q_init)

    return sto.TOProblem(m, c, τ, E, σ_c, σ_t, l, ff, Γ, spDelta, eK, x_init)
end
