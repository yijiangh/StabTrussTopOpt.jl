using StabTrussTopOpt
sto = StabTrussTopOpt

import JSON
using Dates
using Printf
using LinearAlgebra: norm

function parse_TO_problem(ground_fp, load_supp_fp)::sto.TOProblemPar
    ndim, n, m, X, T, B = parse_TO_ground_json(ground_fp)
    F, Bound = parse_TO_support_load_json(load_supp_fp)

    dof_map = zeros(Int, m, 2*ndim) # e_id -> nodal dof
    id_minus = collect(ndim-1:-1:0)
    for i=1:m
        dof_map[i, 1:ndim] = ndim*T[i,1] * ones(ndim) - id_minus
        dof_map[i, ndim+1:2*ndim] = ndim*T[i,2] * ones(ndim) - id_minus
    end

    return sto.TOProblemPar(ndim, n, m, X, T, B, F, Bound, dof_map)
end

function parse_TO_ground_json(file_path::String)
    data = Dict()
    open(file_path, "r") do f
        data_txt = read(f, String)
        data = JSON.parse(data_txt)
    end
    ndim = data["dimension"]
    n = data["node_num"]
    m = data["element_num"]

    X = Matrix{Float64}(undef, n, ndim)
    T = Matrix{Int}(undef, m, 2)
    B = Matrix{Float64}(undef, m, ndim*2)

    # get node coord
    for i=1:n
        X[i,:] = data["node_list"][i]["point"]'
    end

    # get element node ids
    for i=1:m
        T[i,:] = data["element_list"][i]'
    end

    # build compatability matrix
    for i=1:m
        u_id = T[i,1]
        v_id = T[i,2]
        # compute length
        L = norm(X[v_id,:] - X[u_id,:])
        # compute direction cosine
        for k=1:ndim
            B[i, 1:ndim] = (X[v_id,:] - X[u_id,:]) ./ L
            B[i, ndim+1:2*ndim] = -B[i, 1:ndim]
        end
    end

    return ndim, n, m, X, T, B
end

function parse_TO_support_load_json(file_path::String)
    data = Dict()
    open(file_path, "r") do f
        data_txt = read(f, String)
        data = JSON.parse(data_txt)
    end

    ndim = data["dimension"]
    n = data["node_num"]

    n_load_nodes = length(data["point_load_list"])
    @assert(n_load_nodes > 0)
    F = zeros(Float64, n, ndim) # load
    for i=1:n_load_nodes
        load_v = data["point_load_list"][i]["node_id"]
        F[load_v, :] = data["point_load_list"][i]["force"]'
    end

    n_supp_nodes = length(data["support_node_list"])
    @assert(n_supp_nodes > 0)
    Bound = zeros(Int, n, ndim)
    for i=1:n_supp_nodes
        supp_v = data["support_node_list"][i]["node_id"]
        Bound[supp_v, :] = data["support_node_list"][i]["fixities"]'
    end

    # TODO: self-weight
    return F, Bound
end
