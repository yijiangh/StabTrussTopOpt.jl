
function TO_define(par::sto.TOProblemPar)
    ndim = par.ndim
    n = par.n
    m = par.m
    B = par.B
    edof_map = par.dof_map

    ndof = ndim * n
    # free dof
    mask_free = findall(i->i==0, reshape(par.Bound', ndof, 1))
    nfree = length(mask_free)

    ff = reshape(par.F', ndof, 1)[mask_free]

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
    spDelta = spzeros(nfree, nfree)
    for i=1:m
        DDD = spBn[i, mask_free]' * spBn[i, mask_free]
        spDelta[:, i] = DDD[:]
    end
end
